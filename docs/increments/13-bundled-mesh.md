# Increment 13 — one mesh file for ParaView: triangles, constraint edges, feature names

Status: **designed, not started.** The user ruled on the three open questions
on 2026-09-24: U1 (a), U2 (a), U3 yes. See "Ruled by the user" below.

Python only. No C++ change, no new binding. No new runtime dependency. It adds
one optional test dependency (`vtk`, ruling 10).

## Why

The user asked for one file that holds the surface, the constraint edges and
each edge's feature bits, and that ParaView reads correctly. Increment 10 wrote
two PLY files (its ruling 3). That split was designed for QGIS/MDAL, and nobody
had opened the output in ParaView. The measurements below show the PLY path
fails in ParaView in two ways. A later user question ("does river.ply contain
edge type info?") found a third problem: no shipped file says what a feature
bit means.

## What was measured

Everything below was run, not recalled. It used `vtk` 9.7.0 from PyPI
(`uv run --no-project --with vtk python …`). Each candidate file was written by
hand, the way our own writer would write it, not by VTK's writers. Each was
read back with the reader ParaView uses for that format:

| format | reader | ParaView uses it for |
|---|---|---|
| `.ply` | `vtkPLYReader` | `.ply` |
| `.vtp` | `vtkXMLPolyDataReader` | `.vtp` |
| `.vtu` | `vtkXMLUnstructuredGridReader` | `.vtu` |
| legacy `.vtk` | `vtkPDataSetReader` (also `vtkPolyDataReader`, `vtkDataSetReader`) | `.vtk` |

The probe points carry arbitrary 17-digit doubles at UTM scale (x≈430 000,
y≈6 900 000). "Exact" means `np.array_equal` against the written array. The
PLY rows show that a probe that should fail does fail.

**The small bundle: 4 points, 2 triangles, 2 constraint edges, a mask per edge.**

| candidate | points dtype | exact | cells read | cell data |
|---|---|---|---|---|
| PLY, `face` + `edge` in one file, ASCII | float32 | **no** | 2 triangles, **0 lines** | none |
| PLY, `face` + `edge` in one file, binary | float32 | **no** | 2 triangles, **0 lines** | none |
| PLY, `edge` only (increment 10's second file) | float32 | **no** | **0** | none |
| VTP, Float64, ASCII and base64 binary | float64 | yes | 2 lines + 2 triangles | `feature_mask` |
| VTU, Float64, ASCII | float64 | yes | 2 lines + 2 triangles | `feature_mask` |
| legacy `.vtk` 4.2, `POINTS … double`, ASCII | float64 | yes | 2 lines + 2 triangles | `feature_mask` |

So `vtkPLYReader` reads vertices and faces only. It drops `edge` elements with
no warning, whether or not faces are present. It always narrows points to
float32, whatever the file declares. At northing 6.9e6 the float32 step is
2²²⁻²³ = 0.5 m.

**Four more findings decide the details.**

1. **VTK orders cells as lines first, then polygons, whatever order the file
   declares them in.** A VTP with `<Polys>` written before `<Lines>` still reads
   back as cells `[line, line, tri, tri]`. Cell data is matched to cells by
   that index. Cell data written triangles-first is therefore silently given to
   the wrong cells: no error, no warning, error code 0. This was measured with
   distinct values `33 44 11 22` against `11 22 33 44`, in XML versions 0.1
   and 1.0.
2. **A cell array must cover every cell.** A `feature_mask` with entries for
   the lines only makes `vtkXMLPolyDataReader` log an error and return an
   **empty dataset** (0 cells, 0 arrays). So triangles must carry a value.
3. **Strings.** In VTP ASCII, a `String` array is stored as byte codes
   (`114 105 118 101 114 0` for `river`), and VTK's own writer does the same.
   The plain words `river road` did not come back as strings. The reader
   returned an **empty dataset**, so a malformed name loses the whole mesh,
   not just the name. Legacy `.vtk` stores strings as plain text lines. It
   percent-encodes the characters that would break a line: VTK's own writer
   emits `+proj=utm%20+zone=33`. Reading back decodes `%20` and `%25`
   (measured: `a%2541%20b` → `a%41 b`). A raw space in a one-string array also
   survived, but the writer will not rely on that.
4. **`SCALARS` versus `FIELD` in legacy cell data.** `vtkPolyDataReader`, by
   default, returns only the **first** `SCALARS` block in a `CELL_DATA`
   section. ParaView's `vtkPDataSetReader` returns all of them. Arrays placed
   in a `FIELD` block inside `CELL_DATA` are returned by all three readers.

Empty cases: `LINES 0 0`, no `LINES` block, and `FIELD features 0` all load
correctly (1 triangle, cell data intact).

**Cost at DEM scale.** A synthetic 1000 × 1000 grid: 1 000 000 vertices,
1 996 002 triangles, 999 constraint edges, with mm-scale coordinate offsets and
17-digit z. The Python writer used `list(map(repr, a.ravel().tolist()))`, which
is not the shipped PLY writer's per-point generator. Times are from one run on
the development Mac and are indicative only.

| file | size | write (Python) | read (VTK) | points exact |
|---|---|---|---|---|
| legacy `.vtk` ASCII | 90.4 MB | 2.3 s | 1.2 s | yes |
| legacy `.vtk` binary (big-endian) | 63.9 MB | 0.06 s | 0.04 s | yes |
| VTP ASCII | 102.0 MB | 1.7 s | 1.4 s | yes |
| VTP binary (base64 inline, Int64 indices) | 127.8 MB | 0.2 s | 0.7 s | yes |
| PLY binary, shipped writer, faces only | 49.9 MB | 0.04 s | 0.2 s | **no (float32)** |
| PLY ASCII, shipped writer, faces only | 86.4 MB | 3.7 s | 0.3 s | **no (float32)** |

Text costs about 1.4× the size of binary. At a million vertices it adds about
2 s to the write and about 1 s to the read. Both scale linearly: expect about
0.9 GB and 20 s to write at 1e7 vertices. A refined TIN is much smaller than
the DEM grid it comes from, so the realistic case is well below the table.

**Dependency facts.** `vtk` 9.7.0 has wheels for cp312, cp313 and cp314 on
macOS (x86_64, arm64), manylinux (x86_64, aarch64) and Windows. They are
103–140 MB; the Linux x86_64 wheel is 139.7 MB (PyPI JSON). It imports on
3.12.12, 3.13.15 and 3.14.7 here. It requires `matplotlib`. It contains no
GDAL or PDAL reader (`hasattr(vtk, "vtkGDALRasterReader")` is `False`).
`tools/check_prohibited_deps.py` passes with `viewer = ["vtk>=9.3"]` added to
`pyproject.toml`. With `"rasterio"` added beside it, the same run fails and
names it, so the pass is a real result. `pyproject.toml` was restored.

To reproduce: the probe scripts are not kept in the tree. Every row above
comes from writing a file in the format that row names and calling the named
reader's `Update()`.

## Rulings

1. **The format is legacy VTK, version 4.2, `DATASET POLYDATA`, suffix
   `.vtk`.** All three VTK candidates carry float64, lines, triangles and cell
   data, so plain text decided it. It is the only one whose metadata a person
   can read with `head`: the feature names, the CRS and the fingerprint appear
   as words. In VTP ASCII they appear as byte codes, and a writer that gets the
   byte-code form wrong loses the whole dataset (finding 3). Of the three VTK
   formats, legacy `.vtk` is also the smallest in ASCII and the fastest to
   read in ASCII. Its
   writer has the same shape as the PLY writer (a line-oriented header and
   bodies), so the code patterns are shared.
   - Version 4.2 was chosen over 5.1 because it is the older layout and every
     VTK reader accepts it. VTK 9.7 was measured reading it.
   - **VTP is rejected.** It is ParaView's native format and compresses, but
     its strings cannot be read in the file, and it is larger in both
     encodings as measured. If binary size ever becomes the problem, VTP with
     zlib is the upgrade path. That would be a second writer and is not needed
     now.
   - **VTU is rejected.** It reads correctly but adds a per-cell type array and
     the unstructured-grid model to a dataset that is a surface plus polylines.
     PolyData is the exact fit.
   - **A single PLY with `face` and `edge` is rejected.** ParaView drops the
     edges and narrows the points (measured). MDAL's 1D/2D caveat from
     increment 10's ruling 3 also still applies.
   - **MDAL does not read `.vtk`.** This file is for ParaView and not for QGIS.
     What happens to the QGIS path is question U1.

2. **Coordinates are `double` end to end.** The header is `POINTS N double`.
   ASCII values are written with `repr`, as increment 10's ruling 2 already
   does, because `repr` is the shortest text that reads back as the same
   double. Binary values are `>f8`. Both were measured bit-exact through
   `vtkPDataSetReader` at 1e6 vertices. Z comes from the caller as an `(N, 3)`
   array, exactly as in increment 10's ruling 4. `--flat` is unchanged.

3. **ASCII is the default. `--binary` asks for packed records.** The user
   cannot read binary files, and asked for text. Increment 10's ruling 1
   argued for binary on parse time; the measured cost of text is about 1.4× the
   size and about 1 s more read time per million vertices. That is not a strong
   enough reason to give a person a file they cannot open. Binary legacy VTK
   must be big-endian by the format's definition, so the writer uses
   `>`-prefixed dtypes. This is the opposite of the PLY writer's `<`, and a
   test pins it (see Tests).

4. **Cell order is fixed: every `LINES` cell comes before every `POLYGONS`
   cell, in the file and in every cell array.** This matches VTK's internal
   order (finding 1). If the file order and the array order disagree, VTK
   gives the values to the wrong cells in silence. `LINES` is always written,
   as `LINES 0 0` when there are no constraint edges; that form was measured to
   load.

5. **Every cell array covers every cell. Triangles carry 0.** An array covering
   only the lines empties the dataset (finding 2), so leaving triangles out is
   not possible. A triangle's `feature_mask` of 0 and an unclassified
   constraint edge's 0 are told apart by cell type, which the file always
   records. We do not add a separate "is this an edge" array, for two reasons.
   The cell type already says it. And any second array has to be kept
   consistent with the cell order, which is the thing ruling 4 exists to
   protect. Future per-triangle data (slope, DEM error) will follow the same
   rule in the other direction: lines get a fill value.

   **How ParaView colours this dataset.** A colour-by selection applies one
   lookup table to all cells of the chosen array. Colouring by `feature_mask`
   paints every triangle, and every unclassified edge, with the colour for 0,
   and paints each constraint edge by its mask value. Colouring by a
   per-feature array such as `river` (ruling 6) paints river edges as 1 and
   everything else as 0. That is the view a person should use. `feature_mask`
   is written as the `SCALARS` block, so it is the active cell scalar. Two
   things were **not** measured, because nobody here has run ParaView's GUI:
   whether ParaView colours by it when the file is opened, and whether lines
   lying in the triangles' plane draw on top of them. Both are items for the
   manual acceptance step.

6. **The edge vocabulary travels in the file, as names and as a fingerprint.**
   This carries out increment 7's ruling 2 (`07-edge-properties.md`, "How a
   producer and a consumer come to agree which bit means 'river'", mechanism
   2). That ruling says the fingerprint is "a field of whatever serialized
   artifact carries a mesh … exactly as CRS is". Increment 10 does not follow
   it (see ruling 9). The file carries:
   - **Dataset `FIELD FieldData`:**
     - `feature_bits` (`unsigned_int`) and `feature_names` (`string`): the
       vocabulary's `(bit, name)` table, sorted by bit.
     - `feature_vocabulary` (`string`): `EdgeVocabulary.fingerprint()`.
     - `elevation` (`string`): today `none (z=0, --flat)`, the same text as
       increment 10's `FLAT_COMMENT`.
     - `crs` (`string`), only when `--crs` was given.

     The table is the whole vocabulary, not only the bits in use, so a reader
     can tell which bits exist and which are merely unset.
   - **Cell data:**
     - `feature_mask` (`unsigned_int`, the `SCALARS` block): the raw word, for
       scripts.
     - A `FIELD features` block with one `unsigned_char` 0/1 array for each
       property that is **set on at least one edge**, named by the property.
       These arrays are for people: a ParaView user picks `river` from the
       colour menu. Increment 7 keeps names out of C++, and that is not
       touched here: the arrays are built in Python from the vocabulary, and
       C++ still carries only the bits.

     Per-feature arrays are written only for properties that occur. With all
     seven `DEFAULT_VOCABULARY` properties they would add about 28 MB of zeros
     to a million-vertex ASCII file. The consequence: the set of per-feature
     arrays depends on the data. **Scripts must read `feature_mask` and the
     table, not the array names.** The file's schema is the table.
   - **The writer takes the `EdgeVocabulary`, not loose name lists.** It
     computes the table, the fingerprint and the per-feature arrays from that
     one object, so the three cannot disagree. It turns each mask into names
     with `vocabulary.names(mask)`, which raises on a bit the vocabulary does
     not name (increment 7's mechanism 3). So a mask from a different
     vocabulary is refused at write time, not written as an unexplained
     number. `features.py` imports only `hashlib` and `pydantic`, so `io/`
     depending on it adds nothing to what `io/` can reach. `io/` still does not
     import `_core`.
   - **Refusal on read** (the other half of mechanism 2) is not in scope. The
     project has no mesh reader. The fingerprint is written so that the first
     reader can refuse a mismatch.

7. **Strings in the file are ASCII, and no line can be forged.** Every string
   value (`crs`, `elevation`, names, fingerprint) goes through one encoder:
   - control characters are refused, as in the PLY writer;
   - non-ASCII is refused, as in the PLY writer;
   - `%` becomes `%25` and space becomes `%20`, which is VTK's own convention
     (finding 3, measured both ways).

   The title (line 2 of the file) is the fixed text `rasputin mesh`. It never
   holds user input, because a title line has no escaping. The CLI's existing
   step that turns `ValueError` into a `--crs` usage error covers this writer
   too.

8. **CLI: the suffix of `--out` chooses the format.**
   - `--out x.vtk` writes the bundle.
   - `--out x.ply` keeps today's PLY behaviour (subject to U1).
   - Any other suffix is a usage error that names both valid suffixes.
   - `--out-edges` together with a `.vtk` `--out` is a usage error, because the
     edges are already in the file.
   - Encoding becomes one flag pair, `--binary/--ascii`, defaulting to ASCII
     for `.vtk`. Whether PLY's default changes too is question U2.

   There is no `--format` flag, because it could contradict the suffix and then
   one of the two would have to win in silence. A suffix is also what ParaView
   itself uses to pick a reader, so the file is always named for the reader
   that opens it.

9. **Increment 10's PLY edge file breaks increment 7's ruling 2.**
   `river-edges.ply` carries `feature_mask` as a bare `uint`: 0 on the outline,
   1 on the river. Nothing in the header says what bit 0 means, and no
   fingerprint is present. `io/ply.py` and the `mesh` command have never
   carried either. This is a shipped defect. It is not an open design
   question. How it is closed depends on U1:
   - If PLY stays, the edge file gets `comment feature_bit <bit> <name>` lines
     and `comment feature_vocabulary <fingerprint>`, built from the same
     `EdgeVocabulary` through the same `names()` check. MDAL ignores comments,
     so this is for people and future readers, as increment 10's CRS comment
     (its ruling 5) already is.
   - If PLY is dropped, the defect goes with it.

   Either way it is closed in this increment's PR, per
   `docs/increments/README.md` ("fixed in that increment's PR, or it is not
   recorded").

10. **Verification: yes, a test reads the output back with VTK.** It uses an
    optional extra and a separate CI step.
    - Add `viewer = ["vtk>=9.3"]` to `[project.optional-dependencies]`. It is a
      test-only extra, not a runtime dependency and not in `dev`.
    - The VTK suite starts with `pytest.importorskip("vtk")`.
    - `.github/workflows/main.yaml` gets a step on every matrix leg, in the
      same shape as "Tests with the codecs extra": install `.[dev,viewer]` and
      run only that suite with `--no-cov`.

    Reasoning:
    - The defect this increment fixes was only found by opening the output in
      the real reader. A suite that only parses our own bytes would have passed
      the float32 narrowing and the dropped edges, because both happen inside
      the reader.
    - Finding 1 (cell data given to the wrong cells in silence) can only be
      caught end to end by the reader that does the reordering.
    - The cost is one 140 MB download per leg, only in that step.
    - Wheels exist for all three CI Pythons (measured on PyPI).
    - `vtk` is not on `CLAUDE.md` §2's list, and the gate agrees (measured
      above).
    - A leg whose `vtk` wheel disappears fails loudly in that step. It does not
      silently skip, because `importorskip` is only reached after the install
      has succeeded.

    ParaView itself cannot run in CI. Acceptance keeps one manual step.

11. **Numbering: this is increment 13, `13-bundled-mesh.md`.** 12 is reserved
    for the raster adapter (`ROADMAP.md`, MVP gap 1). `10b` is not used,
    because 10 has shipped (#86): this increment reverses one of its rulings
    and adds a format, which is new work and not a sub-step of 10. The branch
    `increment13-bundled-mesh` already matches.

## Ruled by the user

Each of these sets two sound principles against each other. Each changes a
shipped ruling of increment 10, and says which one. **The user chose U1 (a):
keep PLY for QGIS and add the vocabulary; U2 (a): PLY is text by default too;
U3: yes, `vtk` as a test-only extra.** The options are kept so the reasons stay
on record.

**U1 — what happens to the PLY path. Changes increment 10's ruling 3 (two
files), and possibly its ruling 1 (the format is PLY).**

- **(a) Keep PLY as the QGIS output, unchanged in shape, and fix ruling 9's
  vocabulary gap in it.** *Recommended.*
  - The principle this keeps: QGIS reads meshes through MDAL, and MDAL reads
    PLY but not VTK (increment 10's ruling 1). This is a geospatial tool, and
    QGIS is the likeliest non-developer reader.
  - The cost: two writers, and a two-file path that the user's direction
    ("bundling these together is a better choice") moves away from.
  - Note that nobody has ever opened the PLY in QGIS either. Increment 10's
    manual acceptance is still unrecorded. So (a) keeps a path whose one
    reason for existing is still unverified.
- **(b) Drop PLY and `--out-edges`. `.vtk` is the only output.**
  - The principle this keeps: one output, one writer, with no path that
    silently loses data. In ParaView today, the PLY path does both: it drops
    the edges and narrows the coordinates.
  - The cost: QGIS can read nothing the engine writes until a QGIS-readable
    format comes back. Rasputin is 0.2.0.dev0 and nothing outside the repo
    consumes the PLY.
  - Saves ~10 lines against (a).
- **(c) Keep PLY but deprecate it.** A warning goes on stderr, and the defect
  is not fixed. This is not recommended: it keeps both costs and removes
  neither.
- Not on offer: bundling edges into the PLY. The measurements rule it out.

The recommendation is (a), with a condition: when someone opens a PLY in QGIS
(increment 10's acceptance step), record the result in `10-mesh-output.md`. If
QGIS turns out not to need or not to read it, take (b) then.

**U2 — the PLY default encoding. Changes increment 10's ruling 1 (binary by
default).** This only matters if U1 keeps PLY.

- **(a) Flip PLY to ASCII by default. `--binary/--ascii` becomes one flag pair
  for both formats.** *Recommended.*
  - The principle: the user reads the output, and one command should not have
    two defaults. `--ascii` keeps working as the explicit spelling.
  - The cost: the shipped writer's ASCII path is the slowest writer measured
    (3.7 s for 2e6 triangles against 0.04 s). The PLY tests that assume a
    binary default must be updated by `@tester`.
- **(b) Keep PLY binary by default and `.vtk` ASCII by default.** This keeps
  increment 10's reasoning where it was made. The cost: the same command
  defaults to different encodings depending on the suffix.

**U3 — the `vtk` test extra (ruling 10).** This is recommended. It is listed
here because it is the first test dependency heavier than the code it tests,
and the user may prefer a manual-only check. Without it, the pure suite still
pins the bytes, but nothing in CI can detect finding 1 or a reader-side
narrowing.

## Prior art in `legacy/`

```sh
grep -rliE "vtk|\.vtp|\.vtu|paraview|FieldData" legacy/
```

returned no files. The legacy wrote XDMF through `meshio`
(`10-mesh-output.md`, "Prior art"). Nothing is carried across, and no
`@migration-expert` step is needed.

## Files

```
src_python/tin_engine/io/vtk_legacy.py   # new. arrays + EdgeVocabulary -> bytes; no path
src_python/tin_engine/io/__init__.py     # re-export write_vtk
src_python/tin_engine/io/ply.py          # U1(a): vocabulary comments; U2(a): default flip
src_python/tin_engine/cli.py             # suffix dispatch, refusals, --binary/--ascii, docstring
pyproject.toml                           # the `viewer` extra (U3)
.github/workflows/main.yaml              # the VTK read-back step (U3)
project_structure.md                     # io/ gains vtk_legacy.py
docs/increments/10-mesh-output.md        # a note under rulings 1 and 3 pointing here
ROADMAP.md                               # row 13; MVP gap 4 text
```

The module is not named `io/vtk.py`. A test file that does `import vtk` next
to a module named `vtk` is easy to misread.

Signature, for `@tester` to write against:

```
def write_vtk(
    vertices,            # (N, 3) float64; z is the caller's
    *,
    triangles,           # (T, 3) uint32
    edges,               # (E, 2) uint32, may be empty
    edge_masks,          # (E,) uint32
    vocabulary,          # tin_engine.features.EdgeVocabulary
    fields=(),           # Sequence[tuple[str, str]], e.g. ("crs", "EPSG:25833")
    binary=False,
) -> bytes
```

It raises `ValueError` in these cases:
- `vertices` is not `(N, 3)`;
- `edge_masks` does not have one entry per edge;
- a mask carries a bit the vocabulary does not name (from `names()`);
- a field name is not `^[a-z][a-z0-9_]*$`, or collides with a reserved name
  (`feature_bits`, `feature_names`, `feature_vocabulary`);
- a string value is not ASCII or has a control character.

Data flow. It is increment 10's flow with one sink instead of two:

```
mesh + z array            -> (N, 3) vertices
mesh.triangles            -> (T, 3) triangles
cli._constraint_arrays    -> (E, 2) edges + (E,) masks    (unchanged)
DEFAULT_VOCABULARY        -> vocabulary
arrays + vocabulary       -> io.vtk_legacy.write_vtk -> bytes   (pure)
bytes + Path              -> cli writes                          (the only I/O)
```

## Tests for `@tester`

**The invariant-critical suite is the cell order (rulings 4 and 5).** It is the
one mutation round worth paying for here. The mutant to kill: the writer emits
`POLYGONS` before `LINES`, with cell data in that same order. That file is
internally consistent, and VTK misreads it in silence. The pure suite must kill
it through the file-order assertion. The VTK suite must kill it through
semantics: line cells map to masks. Both must be shown red against it.

Pure suite (no VTK; a parser written in the test, as increment 10 did):

- Header: `# vtk DataFile Version 4.2`, title `rasputin mesh`, `ASCII` or
  `BINARY`, `DATASET POLYDATA`, `POINTS N double`.
- Points round-trip bit-exact. Include 430000.001 and 6900000.001, which
  float32 cannot hold, plus 17-digit random values.
- `LINES` comes before `POLYGONS`. `CELL_DATA` counts `E + T`. The first `E`
  entries of `feature_mask` are the masks, and the last `T` are 0.
- One per-feature array for each property set on at least one edge, and none
  for the others. Each value is the bit's state. A mask with no bits set
  produces no `FIELD features` block.
- Dataset field data: bits and names equal the vocabulary sorted by bit, and
  the fingerprint equals `vocabulary.fingerprint()`. `crs` and `elevation`
  appear when given.
- A mask with an unnamed bit raises. A reserved field name raises.
- Strings: space and `%` are encoded. Control characters and non-ASCII are
  refused, naming the character.
- Zero edges: `LINES 0 0` and a well-formed file.
- ASCII and binary decode to the same arrays.
- Binary is big-endian. Use a value whose little-endian reading differs, and
  assert on the raw bytes.

VTK suite (`importorskip("vtk")`, run by ruling 10's CI step):

- Read with `vtkPDataSetReader` and with `vtkPolyDataReader`, in both
  encodings.
- Points are float64 and bit-exact.
- Cell `i < E` is a line (VTK type 3) whose point ids are edge `i`, and whose
  `feature_mask` is mask `i`. Cells `i ≥ E` are triangles (type 5) with 0.
- The per-feature arrays come through both readers. This relies on the `FIELD`
  placement (finding 4).
- The FieldData strings round-trip, including a CRS containing a space and a
  `%`.

CLI suite:

- `.vtk` writes one file.
- `--out-edges` with `.vtk` is refused.
- An unknown suffix is refused.
- `--binary` writes the `BINARY` header.
- `--crs` refusals surface as usage errors for `.vtk` too.
- Existing `.ply` tests stay green, updated only as U1 and U2 require.
- `road-crosses-river` through the CLI, then the VTK reader, gives 9 points,
  12 triangles and the constraint lines. Under U1(a) or U2(a), a PLY edge file
  carries the vocabulary comments.

## LOC estimate

This is an estimate, not a measurement. It is counted in `CLAUDE.md` §2's unit
(comments and docstrings excluded).

| file | what | est. |
|---|---|---|
| `io/vtk_legacy.py` | header, points, lines, polygons, cell and dataset field blocks, ASCII and binary bodies, string encoder, validation | ~105 |
| `io/__init__.py` | re-export | ~2 |
| `cli.py` | suffix dispatch, two refusals, `--binary/--ascii`, vocabulary and fields passed through | ~35 |
| `io/ply.py` | U1(a): vocabulary comments (~10); U2(a): default flip (~2) | ~12 |
| | **total** | **~155** |

Increment 10's estimate overran by 39% (`10-mesh-output.md`,
"Reconciliation"). On the same bias this lands near 215, about a third of the
700 ceiling. The `pyproject.toml` and CI lines (~8) are configuration, not
production code. Under U1(b) the total drops by ~10, and the PLY module and its
suite are deleted: the deletion counts as zero added lines.

## Acceptance

```sh
rasputin mesh road-crosses-river --flat --out /tmp/rr.vtk
head -40 /tmp/rr.vtk          # readable: names, fingerprint, coordinates
```

A person opens `/tmp/rr.vtk` in ParaView and records here the version and:
- that the surface and the constraint lines both show;
- that colouring by `river` isolates the river edges;
- what colouring by `feature_mask` looks like on load (ruling 5's two
  unmeasured items).

## Not in scope

- A mesh reader, and therefore the read-side fingerprint refusal.
- VTP, compression, and time series.
- Any QGIS-readable format other than the existing PLY.
- Per-triangle datasets.
- Any change to `_core`, the bindings, `viz/` or `features.py`.
