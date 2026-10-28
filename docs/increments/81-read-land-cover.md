# Increment 81: `read_mesh` returns the land-cover code of each triangle

**Status: built, approved** (red 03df40e9, green 1f8cd59c, red f88d5cf7, fix bbd5198d, 31 counted lines); designed by `@architect`, 2026-10-20, on `8fd1e118`, master,
branch `worktree-mesh-landcover`. **Short path**
(`docs/increments/README.md@aabac548:76`): about 25 counted lines, no
C++ file, no refine or mesh code.
Number 81 because 80 is the highest number found in `docs/increments/` over
every local branch and remote-tracking ref (`git ls-tree --name-only <ref>
docs/increments/`) and every worktree under `.claude/worktrees/`, checked
2026-10-20.

## 1. What changes

`rasputin mesh --land-cover` writes each triangle's land-cover code into the
`.vtk` as the cell array `land_cover_code` (`int`, every line cell 0, lines
first; `src_python/tin_engine/io/vtk_legacy.py@8fd1e118:164`) and what the
codes are as the dataset string `land_cover_codes`. `read_mesh`, the reader
the core offers plugins (`tin_engine.plugin_api`), stops after `POLYGONS`
(`src_python/tin_engine/io/vtk_read.py@8fd1e118:84`) and so drops both. A
reader of rasputin's own files should return what the writer writes; this
increment makes it return the land-cover codes.

### 1.1 Interface

`MeshArrays` (frozen dataclass, `src_python/tin_engine/io/vtk_read.py`)
gains two fields after `crs`, both with defaults so that
`MeshArrays(vertices, triangles, crs)` stays valid:

| Field | Type | Value |
|---|---|---|
| `land_cover_code` | `npt.NDArray[np.int32] \| None = None` | `(T,)`, one per triangle, in `triangles`' order: the file's `land_cover_code` array with its line-cell entries dropped. `None` when the file has no such array |
| `land_cover_codes` | `str = ""` | the dataset string `land_cover_codes`; `""` when the file has none |

`int32` because that is what the file holds (`int`); no widening.
`read_mesh(data: bytes) -> MeshArrays` keeps its signature.

### 1.2 Behaviour

- After `POLYGONS`, the reader walks the rest of the file by keyword:
  `POINT_DATA n` and `CELL_DATA n` set the current count; `SCALARS name
  type 1` is followed by its `LOOKUP_TABLE` line and a body of count values,
  skipped; `FIELD name k` is followed by `k` arrays (`name 1 count type`),
  each skipped by its size unless it is `land_cover_code`. Sizes come from
  the existing `_SIZE` table; nothing past `CELL_DATA` other than
  `land_cover_code` is kept.
- The dataset string `land_cover_codes` is read in the existing loop over
  the dataset's field data, beside `crs`.
- **Absent array:** `land_cover_code` is `None` and `land_cover_codes` is
  whatever the file holds (`""` if no string). Not an error: a mesh without
  `--land-cover` is a valid mesh.
- **Refusals** (`ValueError`, plain words, beside the existing ones):
  - a `land_cover_code` array whose count is not lines + triangles: "the
    mesh's land_cover_code array has N entries for M cells";
  - a `land_cover_code` array whose type is not `int`: "the mesh's
    land_cover_code array is <type>, not int". Same-size types (`float`,
    `unsigned_int`) would otherwise be read as wrong numbers in silence.
  - a file that ends inside any of the new blocks: the cursor's existing
    "the mesh file is truncated".
- The line cells' entries are dropped by position (lines come first, as the
  writer guarantees); their values are not checked.

The `read_mesh` row of `docs/plugins.md` (the plugin interface table) is
updated in the green commit to name the two new fields.

## 2. Tests `@tester` writes first

In `tests/python/test_io_vtk_read.py`, through its `mesh_bytes` helper
extended with `triangle_codes`, `land_cover_codes` and `class_fractions`.
Fixture, as measured on 2026-10-20 with `write_vtk(..., binary=True)`: the
module's 9 vertices and 5 triangles, edges `[[0, 1], [1, 2]]` with mask 1,
codes `[142, 322, 332, -7, 2**31 - 1]` give `CELL_DATA 7`,
`FIELD features 2` (`river`, then `land_cover_code`) and the header
`land_cover_code 1 7 int` as the file's last array, 53 bytes from header to
end of file.

1. **Round trip with lines:** the codes above and the string
   `"CORINE Land Cover level-3 code"` come back equal, dtype `int32`, shape
   `(5,)`; vertices, triangles and crs bit for bit as today. The negative
   code and `2**31 - 1` catch an unsigned or wrong-width read.
2. **Round trip without lines** (no `LINES` block): the same codes back.
3. **Absent:** a mesh written without `triangle_codes` reads with
   `land_cover_code is None` and `land_cover_codes == ""`.
4. **Codes without a string:** `triangle_codes` given, `land_cover_codes`
   empty: codes back, `land_cover_codes == ""`.
5. **Other arrays around it:** with `class_fractions` (numeric dataset
   arrays) and an edge mask setting two vocabulary bits (two feature arrays
   before `land_cover_code`), the codes still come back right.
6. **Wrong count refused:** the fixture's bytes with the header rewritten to
   `land_cover_code 1 6 int` and one 4-byte value cut: `ValueError` naming
   `land_cover_code`.
7. **Wrong type refused:** the header's `int` rewritten to `float` (same
   size): `ValueError` naming the type.
8. **Truncated:** the fixture's bytes cut inside the `land_cover_code`
   body: `ValueError` "truncated".

No mutation round: the suite is not invariant-critical (a value reader,
`docs/increments/README.md`, *Cost constraints*).

## 3. Left out

- The other arrays the writer puts past the geometry: the per-point
  `elevation` (equal to the vertices' z), the per-cell `feature_mask` and
  per-feature 0/1 arrays, the land-cover fraction table
  (`land_cover_offsets`, `land_cover_classes`, `land_cover_fractions`) and
  the feature vocabulary. Each is its own step when a reader needs it.
- The constraint edges themselves (`LINES`), still skipped.
- ASCII files, still refused.

## 4. Estimate

About 25 counted lines (`CLAUDE.md` §2): two fields, the dataset-string
branch, a keyword walk of about 15 lines, two refusals. Under 50: short path.

## Prior art: legacy and literature

*Literature.* The file layout is VTK's legacy format, version 4.2
(Kitware, *VTK File Formats*, "Simple Legacy Formats": `POINT_DATA`,
`CELL_DATA`, `SCALARS`, `FIELD`; binary data big-endian). The reader follows
it for the keywords rasputin's writer emits and no others. No method, no
novelty claim.

*Legacy.* Nothing to carry: the legacy tree wrote cell data through meshio
and never read a `.vtk` back.

```
$ git grep -n -i -E "land_cover_code|CELL_DATA" legacy-archive -- legacy
legacy-archive:legacy/rasputin/writer.py:29:            writer.write_data(time, cell_data={"triangle": {"shade": shade}})
legacy-archive:legacy/rasputin/writer.py:36:        self._cell_data = dict()
legacy-archive:legacy/rasputin/writer.py:60:    def add_cell_data(self, data, name):
legacy-archive:legacy/rasputin/writer.py:61:        self._cell_data[name] = np.asarray(data)
legacy-archive:legacy/rasputin/writer.py:71:                                   cell_data={"triangle": self._cell_data},
legacy-archive:legacy/rasputin/writer.py:91:                w.add_cell_data(data, name)
```

*Ideas file.* No entry in the private ideas file takes up or overlaps this
step.

## Review

- **Round 1 (short path), 2026-10-20.** Range `bd887652..1f8cd59c` (design, red 03df40e9, green 1f8cd59c), 30 counted lines (estimate about 25), @reviewer: CHANGES REQUESTED. A mesh written with an empty edge vocabulary raised IndexError on read; the plugins.md row for `read_mesh` was not updated. Red-step pins (exact refusal messages, test 5's fixture) and @developer's choices beyond the design accepted.
- **Round 2 (short path), 2026-10-20.** Range `1f8cd59c..bbd5198d`, 31 counted lines in all (estimate about 25), @reviewer: APPROVED. An empty edge vocabulary now reads back (red f88d5cf7, fix bbd5198d); the plugins.md row names the two new fields.
