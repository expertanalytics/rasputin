# Python audit PR D: `audit-encoders` (F7)

Status: **code review round 4 asked for one phrase in this paragraph**
(section Review); code review round 2 approved the code at `083712d`. PR C
(`worktree-audit-geojson`), which this branch was stacked on, merged as
#201. Branch `worktree-audit-encoders` then merged master twice: at
`15f41a82` (master up to `01d98c2b`, #201, which is C) and at `5da8f0f8`
(bringing in #204). The master merge needed one test call fixed,
`4e261dbe` (`tests/python/test_cli_mesh_vtk.py` calls `invoke("mesh", ...)`
as master's shared test driver requires). Net production lines stay +7
(`python3 tools/count_loc.py origin/master HEAD`). Next: a short code review
round 5, then the push on Ola's yes. Follow-up, not blocking: `cli.py` calls
`file_fields(record)` twice on one line. The finding is
`docs/increments/python-audit.md`, F7; its section 6 row D points here.

Citations: into `io/ply.py`, `io/vtk_legacy.py` and `run_record.py`,
pinned to `44fa7f5`; into `cli.py`, `crs.py`, `features.py` and
`fetch/nve.py`, which PR B or C changed, pinned to `b26beb8`. Both commits
are on master (`git merge-base --is-ancestor <commit> origin/master`),
since C was merged rather than rebased, so the pins resolve.

## 1. What F7 said, re-measured

F7 listed five encoders in `cli.py` to move into `io/`, two copies of the
mesh writers' text gate and of their int32 check, two copies of the sorted
vocabulary table, and two copies of the ASCII escape; about -30 net lines.
Read against `b26beb8`:

| F7's item | At `b26beb8` | In PR D |
|---|---|---|
| the fetched station files | not encoded in `cli.py`: `fetch_station_set` returns each file's bytes (`src_python/tin_engine/fetch/nve.py@b26beb8:125-128`); `cli.py` only writes them (`src_python/tin_engine/cli.py@b26beb8:1334-1335`) | nothing to move; F7 was wrong here |
| `palette`'s JSON | `json.dumps(paraview_preset(table, title), indent=1) + "\n"` (`src_python/tin_engine/cli.py@b26beb8:1101`); the document is built in `palettes.py` | stays (section 4, "What stays in `cli.py`") |
| `summary.json` | `json.dumps(dump, indent=2) + "\n"` of a model dump (`src_python/tin_engine/cli.py@b26beb8:2180`) | stays |
| `results.csv` | `csv.writer(...).writerow(...)` of `StationResult`'s fields, through `_cell` (`src_python/tin_engine/cli.py@b26beb8:1984`, `:1991`, `:2026-2030`) | stays |
| the `.ply` comment lines | `_comments(record)` (`src_python/tin_engine/cli.py@b26beb8:1061-1063`) | moves to `io/ply.py` as `field_comments` |
| the text gate | `src_python/tin_engine/io/ply.py@44fa7f5:114-131` and `src_python/tin_engine/io/vtk_legacy.py@44fa7f5:192-200`, two wordings (below) | one copy, `io/mesh_checks.py` |
| the int32 check | `src_python/tin_engine/io/ply.py@44fa7f5:176-177`, `src_python/tin_engine/io/vtk_legacy.py@44fa7f5:114-115` | one copy, `io/mesh_checks.py` |
| the sorted vocabulary table | `src_python/tin_engine/io/ply.py@44fa7f5:108`, `src_python/tin_engine/io/vtk_legacy.py@44fa7f5:108`, and a third copy F7 missed, `fingerprint`'s (`src_python/tin_engine/features.py@b26beb8:163`) | `EdgeVocabulary.table()` |
| the ASCII escape | `cli._ascii` (`src_python/tin_engine/cli.py@b26beb8:1471-1473`), `crs_label` (`src_python/tin_engine/crs.py@b26beb8:112`) | `run_record.escaped_ascii` |

**The net is about +7, not -30.** A prototype of this design, in a scratch
clone of `b26beb8` (not committed), counts +7
(`python3 tools/count_loc.py b26beb8 <prototype>`; section 7). The
duplicated checks are six or fewer code lines a copy, and the shared module
costs its imports. The -30 assumed the encoders move out of `cli.py`;
measured, moving them costs lines (section 4). What PR D buys is one copy
of each rule a file format depends on, so the GeoPackage writer (section 5)
calls them rather than copying them a third time. Question 1.

## 2. Prior art: legacy and literature

*Literature.* No method; three format specifications, each already followed
by the code PR D touches, and nothing new claimed. PLY (Turk, 1994): a
line-oriented ASCII header with `comment` lines, so a line terminator in a
comment forges a header line, as `io/ply.py`'s guard records
(`src_python/tin_engine/io/ply.py@44fa7f5:115-119`). Legacy VTK (Kitware,
*VTK File Formats*): string arrays in `FIELD` data, which `vtk_legacy.py`
escapes (`%` then space) and gates (its ruling 7,
`src_python/tin_engine/io/vtk_legacy.py@44fa7f5:29-30`). GeoPackage (OGC 12-128, Table 1,
"GeoPackage Data Types"): `TEXT` is stored in the database encoding (UTF-8
or UTF-16) and `MEDIUMINT` is a 32-bit signed integer; so a GeoPackage
writer needs the int32 check if it declares the land-cover code
`MEDIUMINT`, and does not need the ASCII gate (section 5).

*Legacy.* `git grep -lE 'backslashreplace|isascii|control char|write_ply|\.ply|DataFile Version|vtk' legacy-archive -- legacy`
returns nothing. The wider
`git grep -liE 'ply|vtk|write_|ascii' legacy-archive -- legacy` returns
`legacy/rasputin/geometry.py`, `legacy/rasputin/mesh.py`,
`legacy/rasputin/reader.py`, `legacy/rasputin/web_visualize.py` and
`legacy/rasputin/writer.py`: the writers there are `meshio` calls
(`write_points_cells`, `XdmfTimeSeriesWriter`, in `writer.py` and
`Mesh.write` in `mesh.py`), with no header text, no gate and no code
check. Nothing to carry across.

## 3. The design

### `io/mesh_checks.py` (new, layer L2, imports nothing first-party)

```python
def checked_ascii(value: str, what: str) -> bytes:
    """`value` as ASCII bytes; a control character (below U+0020, or U+007F)
    or a non-ASCII character is a ValueError naming `what` and the character."""

def check_int32(values: npt.NDArray[Any], name: str) -> None:
    """ValueError `<name> must fit int32` unless every value is in
    [-2**31, 2**31 - 1]; an empty array passes."""
```

The refusals, one wording for both writers:

- `a <what> may not contain a control character; got <repr>`
- `a <what> must be ASCII; got <repr> in <repr of value>`
- `<name> must fit int32`

`io/ply.py` calls `checked_ascii(comment, "comment")` for each comment and
`check_int32(values, "face_codes")`; `io/vtk_legacy.py` calls
`checked_ascii(value, "string")` in `_strings` (its `_checked` goes) and
`check_int32(codes, "triangle_codes")`. Each keeps its own length checks,
whose wordings differ and stay.

A module of its own, not `features.py` (which both writers already
import): `features.py` is the edge-property vocabulary and the feature
shapes that carry it (it imports `hashlib`, `typing`, `pydantic` and
`shapely.geometry`, `src_python/tin_engine/features.py@b26beb8:58-62`),
and a rule for the text a file header may hold is not about edges. Not
`io/__init__.py`: it imports `ply` and `vtk_legacy`, so their importing it
back is a cycle.

### `io/ply.py`

```python
def field_comments(fields: Iterable[tuple[str, str]]) -> list[str]:
    """A file's `(name, value)` fields as header comments, `name value`."""
```

`cli._comments(record)` becomes `field_comments(file_fields(record))`.
`write_ply`'s signature does not change (its `comments` stays), so the
tests that call it keep their calls.

### `features.EdgeVocabulary.table()`

```python
def table(self) -> list[tuple[int, str]]:
    """The (bit, name) pairs in ascending bit order, as files write them."""
```

`fingerprint`, `write_ply` and `write_vtk` call it. The order is the one
`names()` already documents; the fingerprint's bytes do not change.

### `run_record.escaped_ascii`

```python
def escaped_ascii(text: str) -> str:
    """`text` with non-ASCII escaped (`\xe9`), as a file field holds it."""
```

`run_record.py` is layer L0 and pure; `Entry.value` already says its text
is ASCII (`src_python/tin_engine/run_record.py@44fa7f5:88`). `cli.py`'s
seven `_ascii` calls call it, and `crs_label` returns
`escaped_ascii(parsed.to_string())`. `crs` (L1) importing `run_record`
(L0) goes down a layer. The escape is the other half of the gate: what a
run puts in a file field is escaped to ASCII first, so the writers' ASCII
refusal is met only by `mesh --crs`, whose text is not escaped (ruling 5:
free text, refused rather than altered).

### `cli.py`

- imports `escaped_ascii`, `field_comments` and `checked_ascii`; `_ascii`
  and `_comments` go.
- The gallery path's `--crs` check
  (`src_python/tin_engine/cli.py@b26beb8:985-992`) stops encoding an empty
  PLY to see whether the writer refuses; it calls the gate on the same
  lines, `for comment in field_comments(file_fields(record)):
  checked_ascii(comment, "comment")`, inside the same `try`. The refusal
  is the same text, as the PLY writer gives it. The comment above the
  `try` says "the gate's refusals" rather than "the writer's".
- The suffix dispatch (`MESH_SUFFIXES`, `if out.suffix == ".vtk" ... else`,
  `src_python/tin_engine/cli.py@b26beb8:1014-1047`) stays, its `_comments`
  call replaced (one line).
- `values["features"] = escaped_ascii(...)` and `values["features_crs"] =
  ...` as two assignments (`src_python/tin_engine/cli.py@b26beb8:955`):
  with the longer name the one-line `values |= {...}` passes 100
  characters and the formatter spreads it over four lines.

## 4. What stays in `cli.py`, and what was rejected

**The palette JSON, `summary.json` and `results.csv` stay.** Each is one
standard-library call (`json.dumps`, `csv.writer`) on a document or a row
built elsewhere (`palettes.paraview_preset`, `BatchSummary.model_dump`,
`StationResult`). The project's own choices in them are `_cell`'s (a tuple
joined by `;`, None empty) and the header's `class` for `station_class`,
which are the row's wording, owned by `StationResult`'s module
(`catchment_batch.py`, layer L4), which PR F (`worktree-audit-catchment`)
is rewriting; moving them now would conflict with F for no lines saved.
A function in `io/` per call would be wider than the call. `json.dumps`
escapes non-ASCII by default, so the palette's locale-encoded
`write_text` writes ASCII whatever the locale. `project_structure.md`'s
`io/` line, which lists these as known exceptions "to move here in a
follow-up", is rewritten to state the rule instead: a document built
outside `cli.py` and serialised there by one standard-library call is
printing, not an encoder. Question 2.

**Rejected: a suffix-to-encoder table in `io/` (`io/mesh_files.py`, with
`MESH_SUFFIXES` and `mesh_encoders(suffix, mesh, fields, *, binary, codes,
codes_text)`).** Prototyped: +46 net for the PR, +39 more than this design,
because the formatter spreads each writer's call over a line per argument
and the module adds two signatures and its imports. It would make the
GeoPackage writer a table row; the mesh run's move out of `cli.py` (PR H,
`audit-mesh-run`, F8) is where that seam belongs, as the run and its file
writing leave `cli.py` together.

## 5. Room for the GeoPackage writer

The ROADMAP row "GeoPackage output" (`rasputin mesh --out X.gpkg`, on
`origin/master`; built after PR D) says the new writer "joins that shared
check instead of a third copy". Without its design being made here, what it
gets from PR D, and what it already has:

- **`check_int32`**, if it declares `land_cover_code` as `MEDIUMINT` (OGC
  12-128, Table 1); with `INTEGER` (64-bit) it needs none.
- **`EdgeVocabulary.table()`**, for a layer of constraint lines carrying
  feature names (the row's open question), in the order the `.vtk` and
  `.ply` files use.
- **Already-escaped field values**: the run record's file fields
  (`file_fields`) are ASCII by `escaped_ascii`, so a GeoPackage table of
  them holds the same text as the other two formats.
- **Not `checked_ascii`**: a GeoPackage `TEXT` is UTF-8 and has no header
  lines to forge, so the gate is not that format's rule. The one free-text
  field, `--crs` on the gallery path, is checked in `cli.py` before any
  writer runs, so a GeoPackage written from the gallery would still refuse
  `--crs` with a control character, as today's two formats do.
- **The suffix**: one entry added to `MESH_SUFFIXES` and one branch beside
  the `.vtk` one in `cli.py`'s encoder list (or a row of PR H's seam, if H
  lands first). The `--out-edges` refusal beside `.vtk`
  (`src_python/tin_engine/cli.py@b26beb8:822-825`) is the pattern for a
  format that holds the edges itself.
- **Bytes, not a path**: the encoder slot is `Callable[[], bytes]`, and
  `io/` writes no file (`src_python/tin_engine/io/__init__.py@b26beb8:9-15`).
  Python's `sqlite3.Connection.serialize()` (3.11 and later) returns an
  in-memory database's bytes, so a writer can build the GeoPackage in
  `:memory:` and keep both. Checked here: on the project's Python 3.14.7
  with SQLite 3.53.4, `serialize()` of a one-table `:memory:` database
  returns 8 192 bytes. Whether a 678 MB file (the Mur, in the row) should
  be held in memory whole is the writer's question.

## 6. Layering table (`tests/python/test_layering.py`)

| Row | At `b26beb8` | After |
|---|---|---|
| `io.mesh_checks` (new, L2) | none | `""` |
| `io.ply` | `features` | `features io.mesh_checks` |
| `io.vtk_legacy` | `features` | `features io.mesh_checks` |
| `crs` | `""` | `run_record` |
| `cli` | as at `b26beb8` | adds `io.mesh_checks` |

The whole Python suite on the prototype's scratch copy (unchanged tests):
5 failed, 5504 passed, 27 skipped; the five failures are exactly these
rows, `test_the_table_names_every_module` (the new module) and
`test_each_module_imports_exactly_its_row` for `cli`, `crs`, `io.ply` and
`io.vtk_legacy`. `run_record` stays `""`. Every new edge goes down (L1 to L0) or sideways
(L2 to L2, L5 to L2). The old firewall rule for the two writers (no
`_core`, one first-party import) is kept in spirit: their imports are
still pure modules, now two.

## 7. Net production lines: about +7

The prototype's count (`count_loc.py b26beb8 <prototype>`), with the
formatter run; `@developer`'s code may differ by a line or two.

| File | Net | What |
|---|---|---|
| `io/mesh_checks.py` (new) | +15 | 3 import lines, `checked_ascii` 8, `check_int32` 3, `__all__` 1 |
| `io/vtk_legacy.py` | -8 | `_checked` (8) goes; the int32 check 2 to 1; `table`'s line and `_strings`' call rewritten; one import (4 added, 12 removed) |
| `io/ply.py` | -4 | the gate 6 to 1; the int32 check 2 to 1; the `table` line goes, its use rewritten; `field_comments` 2; one import (6 added, 10 removed) |
| `features.py` | +1 | `table` +2, `fingerprint` -1 |
| `run_record.py` | +2 | `escaped_ascii` |
| `crs.py` | +1 | one import |
| `cli.py` | 0 | `_ascii` 2 and `_comments` 2 go; the `--crs` check 1 to 2; the `features` line 1 to 2; `fields, comments = ...` and seven `_ascii` calls rewritten; two import lines (14 added, 14 removed) |

No packing under `# fmt: skip` is designed.

## 8. What the bytes do: the probe

`docs/increments/python-audit-probes/encoder_bytes.py` (its docstring says
how to run it on a scratch copy) prints one line per case: the sha256 of
every file a writer or `rasputin mesh` writes, or the refusal. Part 1 calls
`write_ply` and `write_vtk` directly (5 meshes: empty, one triangle, two
with constraint edges, a NaN z, a UTM easting; text and binary; 8 comment
and field texts: plain, non-ASCII, `\r`, tab, DEL, NUL, `%` and space,
empty; 7 code arrays: none, small, at the int32 bounds, one past each
bound, NaN in a float array, one short; an unnamed mask bit; a reserved
field name) and `crs_label` (an EPSG code, `OGC:CRS84`, a PROJ string,
ASCII and non-ASCII WKT). Part 2 runs `rasputin mesh` in a child process:
the `river` gallery fixture with `--flat` to `.vtk` and to `.ply` with
`--out-edges`, text and binary, with `--crs` plain, non-ASCII, with `\r`
and empty; the two suffix refusals; and a DEM run (a hand-written 21 x 17
GeoTIFF named `dém.tif`, CORINE-coded `--features` named `lånd.geojson`)
to `.vtk`, to `.ply` with `--out-edges`, and to `.ply` with `--record`,
whose escaped names it prints.

Its output at the base is committed beside it,
`encoder_bytes-b26beb8.txt` (124 lines, no `UNCAUGHT`; two runs at the
base are identical). **The expected diff at PR D's head is exactly the 12
lines that carry the PLY control-character refusal**, which reads
`a comment may not contain a control character` instead of `... control
characters`: 8 from `write_ply` (`\r`, tab, DEL, NUL; text and binary) and
4 from `rasputin mesh --crs` with `\r` (`.vtk` and `.ply`, text and
binary). Every hash, every other refusal, the `crs_label` texts and the
`--record` names are unchanged. The prototype's output differs from the
base's in exactly those 12 lines. The probe can fail: a planted double
space in `field_comments`'s output changed 12 other lines (every `.ply`
hash from `rasputin mesh` and the non-ASCII `--crs` refusals).

The green commit's author runs the probe at the head and commits its
output as `encoder_bytes-<head>.txt`; `diff` of the two files is then the
record of what PR D changed.

## 9. Red tests (`@tester`, one commit, before any code)

Lean: no mutation round.

1. **`tests/python/test_io_mesh_checks.py`** (new; red: the module does not
   exist): `checked_ascii("crs EPSG:25833", "comment") ==
   b"crs EPSG:25833"`; `checked_ascii("", "string") == b""`; each of
   `\x00`, `\t`, `\n`, `\r`, `\x1b`, `\x7f` refused with the wording above,
   its `what` and the character's `repr`; a degree sign refused with
   `must be ASCII`, the character and the value. `check_int32`:
   `[2**31 - 1, -(2**31)]` and an empty array pass; `[2**31]` and
   `[-(2**31) - 1]` refused with `<name> must fit int32`.
2. **`tests/python/test_io_ply.py`**: the control-character refusal reads
   `a comment may not contain a control character; got '\r'` (red: the
   plural today); `field_comments([("crs", "EPSG:25833"), ("heights",
   "x")]) == ["crs EPSG:25833", "heights x"]` (red).
3. **`tests/python/test_features.py`**: `DEFAULT_VOCABULARY.table()` is
   the sorted `(bit, name)` pairs (red), and `fingerprint()` equals its
   value at `b26beb8`, written as a literal (green at red: a guard, not a
   red test).
4. **`tests/python/test_run_record.py`**: `escaped_ascii("dém; ø") ==
   "d\\xe9m; \\xf8"` (red). Append at the end of the file:
   `docs/increments/27-node-sampling.md` line 315 cites its line 443
   unpinned.
5. **`tests/python/test_cli_mesh_vtk.py`**, beside
   `test_a_bad_crs_is_a_usage_error_and_no_file`: `--crs "EPSG:25833\rforged"`
   exits 2 naming `--crs` with `a control character` (red).
6. **`tests/python/test_layering.py`**: the table above (red until the
   imports move).

The rest of the suite stays green unchanged: the existing `match="control
character"` and `match="int32"` assertions match both wordings.

### Pinned by the red step (`8dcfa2a`), ruled

`@tester` pinned ten behaviours the design did not word; `@architect` rules
on each. "As today" was checked against both writers at `b26beb8`
(`src_python/tin_engine/io/ply.py@b26beb8:114-131`,
`src_python/tin_engine/io/vtk_legacy.py@b26beb8:192-200`; the int32 checks
at `src_python/tin_engine/io/ply.py@b26beb8:176-177` and
`src_python/tin_engine/io/vtk_legacy.py@b26beb8:114-115`). All ten are
confirmed; none is overruled.

1. **Confirmed.** A value with both a control character and a non-ASCII
   one is refused as a control character: both writers test for control
   characters first today.
2. **Confirmed, read as "the first of the kind refused".** Both writers
   take `next(...)` over the value, so the control refusal names the first
   control character and the ASCII refusal the first non-ASCII one; in
   `"\N{DEGREE SIGN}\r"` the refusal names `\r`, not the earlier degree
   sign. The tests assert exactly that.
3. **Confirmed.** A C1 control (U+0080 to U+009F; U+0085 tested) is not
   below the space nor DEL, so today's control test passes it and the
   ASCII test refuses it, as `must be ASCII`.
4. **Confirmed.** Space and `~` pass; they bound the accepted range.
5. **Confirmed.** NUL, tab, `\n`, `\r`, ESC, U+001F and DEL, each with
   `what` "comment" and "string": the design's six plus U+001F, the last
   character below the space.
6. **Confirmed.** `check_int32` judges an unsigned array by its values:
   uint64 `2**31 - 1` passes and `2**31` is refused. Today's expression
   does the same (run here on the worktree's venv: uint64 `[0, 2**31 - 1]`
   passes, `[2**31]` is refused). The refusal names the `name` given.
7. **Confirmed.** `field_comments` takes any iterable of pairs, a one-shot
   iterator included (the design's signature is `Iterable`); no fields
   gives no comments.
8. **Confirmed.** `EdgeVocabulary.table()` is sorted by bit whatever the
   declaration order (all three copies sort today), and an empty
   vocabulary gives `[]`.
9. **Confirmed.** `escaped_ascii` escapes characters beyond Latin-1 too
   (`\u2603`, `\U0001f600`), as `backslashreplace` does in both copies
   today, and leaves ASCII alone, control characters included: refusing
   those is the gate's job.
10. **Confirmed.** The `--crs` test runs for `.vtk` and for `.ply` and
    asserts `a comment may not contain a control character; got '\r'`,
    exit 2, no file. Today the gallery path's check is a `.ply` write
    whatever the suffix (`src_python/tin_engine/cli.py@b26beb8:985-992`),
    so the `.vtk` run already gets the comment wording; PR D changes only
    the plural, as section 8's 12 lines say.

**`test_layering.py`: 7 red rows, not 5, as intended.** Section 6's five
were counted on the prototype, where `io/mesh_checks.py` existed and the
tests were unchanged. At the red commit the table has the new row and the
module does not exist, so two more fail for that reason alone:
`test_each_module_imports_exactly_its_row[io.mesh_checks]` and
`test_no_module_imports_by_name`. Both go green when the module exists;
`@tester`'s scratch copy with the design applied passed the whole suite.

**Known gap, kept: NaN in a float code array.** `check_int32`, like both
writers today, compares with `min()` and `max()`, and every comparison
with NaN is false, so `[nan, 1.0]` passes and the writer casts the NaN to
an integer without a word. Section 8's probe shows it at the base: the
four "float nan" cases are file hashes, not refusals
(`docs/increments/python-audit-probes/encoder_bytes-b26beb8.txt@8dcfa2a:83-86`).
On this Mac the cast gives 0 (NumPy's `RuntimeWarning: invalid value
encountered in cast`, then 0, both for `astype(np.int64)` and for a `<i4`
field); another platform's cast may give another integer. **Ruling: PR D
keeps it.** PR D is a refactor whose probe must change only the 12 wording
lines; refusing NaN would change four more. No caller in the run passes a
float code array (land-cover codes are integers from the raster), so it
meets only a direct caller of the writers. Default for the later fix: in
its own increment, `check_int32` refuses any value that is not a finite
whole number with the same `<name> must fit int32`, red test first, and
the probe's four lines become refusals.

## 10. `@perf`

Nothing owed. No refine or mesh code changes; the writers run after the
mesh is made, and the probe shows their bytes unchanged. The gate now runs
the same per-character test it ran before, once per comment or string (a
few dozen per file); `check_int32` is the same two reductions.
`tools/bench.py` patches `refine` in `cli`'s namespace, which PR D does not
touch.

## 11. Citations this PR moves, pinned now

Found by running `python3 tools/check_citations.py --base b26beb8` in the
prototype clone, which edits the seven production files (14 at risk).
Pinned in this commit to `44fa7f5`, each quotation re-read there:
`docs/increments/25-plain-output.md` lines 88 (`io/vtk_legacy.py:123-126`,
the three vocabulary arrays), 96 and 273 (`io/vtk_legacy.py:118-119`, the
`land_cover_codes` field), 271 (`io/vtk_legacy.py:124-125` and
`io/ply.py:111`, `feature_bits`, `feature_names` and the `feature_bit`
comments) and 272 (`io/vtk_legacy.py:126` and `io/ply.py:112`, the
fingerprint). Left alone, as PRs B and C left them: unpinned citations
into `cli.py` whose text had already moved (`15f-edge-strip.md` lines 1539
and 1904, `24-release-hardening.md` lines 515 and 517,
`27-node-sampling.md` line 400, `29-nve-reference-catchments.md` lines
3098 and 3114). No unpinned citation points into
`project_structure.md`, `test_layering.py`, `test_io_ply.py` or
`test_features.py`; `test_run_record.py`'s is red test 4's.

`project_structure.md`'s rows for `io/` (the exceptions line),
`io/ply.py`, `io/vtk_legacy.py`, `features.py`, `crs.py` and
`run_record.py`, and a new `io/mesh_checks.py` row, are rewritten by
`@architect` after the green commit, so they describe the code as written.

## 12. Risks

- **Merge order.** PR C lands first (this branch is on it). PR F
  (`audit-catchment-shared`) also edits `cli.py` and `test_layering.py`
  (`git diff --stat b26beb8...worktree-audit-catchment` names both), and
  PR A (`audit-lattice`, being designed) may; whichever of them and D merges second
  runs `git merge-tree --write-tree <its head> master` and resolves every
  file it names, as section 10 of `python-audit.md` rules for C and F.
  The `crs` row and the `cli` row of the table are the likely hunks.
- **The one wording change** reaches a user only through `mesh --crs`
  with a control character on the gallery path; the DEM path does not
  take `--crs` (`--dem records the DEM's own CRS`).
- **`crs` now imports `run_record`.** `run_record` imports only `json`
  and `dataclasses`, so no cycle; `test_layering.py`'s check 3 confirms
  the edge goes down.

## 13. Questions for Ola (defaults hold until he answers)

1. **PR D adds about 7 lines rather than removing 30.** It puts the two
   mesh writers' checks (text that may go in a file header, codes that
   must fit 32 bits) and the ASCII escape in one place each. The main
   reasons: the two copies have already drifted (the `.ply` and `.vtk`
   writers word the control-character refusal differently today), and
   the `mesh --crs` check stops writing an empty `.ply` file in memory
   just to see whether the writer refuses. A lesser reason: the
   GeoPackage writer can use the shared int32 check and the vocabulary
   table instead of copying them (section 5; it needs nothing else from
   PR D). Go ahead at +7, or drop PR D and leave the copies where they
   are? Default: go ahead.
2. **The palette file, `summary.json` and `results.csv` stay written in
   `cli.py`**, each one standard-library call on data built elsewhere;
   moving them would add lines and, for `results.csv`, collide with PR F.
   Default: they stay, and `project_structure.md` says why. If they stay,
   the ROADMAP's "GeoPackage output" row (on `origin/master`), which says
   PR D moves the file encoders out of `cli.py` into `io/`, is wrong about
   PR D and is corrected when D lands (not edited here).
3. **One refusal changes by one word:** `rasputin mesh river --flat --crs
   "x<carriage return>y" --out a.ply` says `a comment may not contain a
   control character` instead of `... control characters`, matching the
   `.vtk` writer's wording. Default: accept.

## 14. As built (`7346a0e`)

**Net production lines: +7** (46 added, 39 removed), measured with
`python3 tools/count_loc.py b26beb8 7346a0e`; every file matches section
7's table: `io/mesh_checks.py` +15, `io/vtk_legacy.py` -8 (4 added, 12
removed), `io/ply.py` -4 (6 added, 10 removed), `features.py` +1,
`run_record.py` +2, `crs.py` +1, `cli.py` 0 (14 added, 14 removed). No
packing under `# fmt: skip`.

**One departure from the design, a comment.** The two long comments in
`io/ply.py`'s old check (`src_python/tin_engine/io/ply.py@b26beb8:115-119`
and `:123-128`) are cut to two lines (`src_python/tin_engine/io/ply.py@7346a0e:115-116`);
its reasons (any control character forges a line, `\r` for CRLF-tolerant
readers, non-ASCII a documented `ValueError` because `--crs` is free text)
are now `io/mesh_checks.py`'s docstring
(`src_python/tin_engine/io/mesh_checks.py@7346a0e:5-9`), where the check
is. Its one historical sentence (the `--crs "x\rcomment forged"` run that
once exited 0) is not carried over; red test 5 now guards that case.
Comments are not counted, so the net is unchanged.

**The probe.** Run at `7346a0e` on a scratch copy (`tools/scratch_copy.py`,
with `PYTHONPATH` taken from its stdout line), its output is
`docs/increments/python-audit-probes/encoder_bytes-7346a0e.txt` (124
lines, no `UNCAUGHT`). `diff encoder_bytes-b26beb8.txt
encoder_bytes-7346a0e.txt` is exactly section 8's 12 lines (39 to 53 odd,
103-104, 115-116), each `control characters` becoming `a control
character`. The base, re-run the same way, reproduces
`encoder_bytes-b26beb8.txt` byte for byte (`cmp`).

**Known gap kept:** NaN in a float code array still passes `check_int32`
(section 9, "Known gap, kept"); the probe's four "float nan" lines are
still file hashes, unchanged from the base.

`project_structure.md`'s rows for `io/` (the exceptions line now states
section 4's rule), `io/ply.py`, `io/vtk_legacy.py`, `features.py`,
`crs.py` and `cli.py` are rewritten against `7346a0e`, and rows for
`io/mesh_checks.py` and `run_record.py` (which had none) are added.

## Review

**PR D (`audit-encoders`), design review, round 1, 2026-10-06.** Range `b26beb8..7fcb47c`. Verdict: CHANGES REQUESTED, prose only: section 3's claim that `features.py` imports only `hashlib` and `pydantic` is false at b26beb8 (`typing` was already there; PR C added `shapely.geometry`, `src_python/tin_engine/features.py@b26beb8:58-62`; section 3 at `docs/increments/python-audit-pr-d.md@7fcb47c:99-103`); section 6 row D's 'Gates beyond review: none' (`docs/increments/python-audit.md@7fcb47c:510`) should list the red tests and the wording change. +7 recounted per file; probe reproduces byte for byte and three planted mutants changed it.

**PR D (`audit-encoders`), design review, round 2, 2026-10-06.** Range `7fcb47c..da99916`. Verdict: APPROVED, prose only: section 3's imports match `src_python/tin_engine/features.py@b26beb8:58-62` (`typing` predates PR C, `3b739ac` added only `shapely.geometry`); row D's gates (`docs/increments/python-audit.md@da99916:510`) match section 9 (`docs/increments/python-audit-pr-d.md@da99916:296-327`); Question 1's drift and in-memory `.ply` reasons hold (`src_python/tin_engine/io/ply.py@da99916:122`, `src_python/tin_engine/io/vtk_legacy.py@da99916:196`, `src_python/tin_engine/cli.py@da99916:990`); `check_citations.py --base b26beb8` exits 0.

**PR D (`audit-encoders`), code review, round 1, 2026-10-06.** Range `da99916..83354aa`. Verdict: CHANGES REQUESTED: `tests/python/test_io_mesh_checks.py@83354aa:10` said 'Committed red: the module does not exist', false since `7346a0e` (red-step scaffolding); line 16's '(the handback lists each)' should point to `docs/increments/python-audit-pr-d.md@83354aa:331-339`. +7 net (46/39), as section 7. Red real at `8dcfa2a` (21 failed, 32 errors, stated reasons), green at `7346a0e`. Suite 5439 passed, 30 skipped on a scratch copy; gates clean. Probe reproduces both committed outputs byte for byte, its diff is the 12 wording lines, and a planted gate mutant changed 12 lines. No `@perf` owed.

**PR D (`audit-encoders`), code review, round 2, 2026-10-06.** Range `83354aa..083712d`. Verdict: APPROVED. Docstrings only: red-step wording past tense at `tests/python/test_io_mesh_checks.py@083712d:10-17` and `tests/python/test_io_ply.py@083712d:495-497`, pointing to `docs/increments/python-audit-pr-d.md@083712d:331`; +7 net (`count_loc.py b26beb8 083712d`); CI to be checked after the push. Non-blocking: `src_python/tin_engine/cli.py@083712d:1018` calls `file_fields(record)` twice.

**PR D (`worktree-audit-encoders`, the shared mesh-writer checks), code review, round 3, 2026-10-07.** Range `4725f12e..4e261dbe` (master merges `15f41a82` and `5da8f0f8`, the `@tester` test fix `4e261dbe`). Verdict: CHANGES REQUESTED. The Status paragraph at `/Users/skavhaug/projects/rasputin/.claude/worktrees/audit-encoders/docs/increments/python-audit-pr-d.md@4e261dbe:3-20` says PR C (`worktree-audit-geojson`) is not pushed and that D gets rebased, but C merged as #201 and D merged master. The code, the merges, the suite and the gates are clean.

**PR D (`worktree-audit-encoders`, the shared mesh-writer checks), code review, round 4, 2026-10-07.** Range `4e261dbe..830eab37`, docs only, +7 net unchanged. Verdict: CHANGES REQUESTED: `/Users/skavhaug/projects/rasputin/.claude/worktrees/audit-encoders/docs/increments/python-audit-pr-d.md@830eab37:7` says merge `15f41a82` brought in "#200 to #203", but it brought in master from #187 up to #201 (`git log --first-parent b26beb83..01d98c2b`). The other claims check out: #201 is C, the pins are on master, `4e261dbe` is the test fix.
