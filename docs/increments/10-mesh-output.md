# Increment 10 — mesh output

Status: design. No production code and no test on this commit.

Fills MVP gap 4 in `ROADMAP.md`: *"Mesh output. There is no writer of any
kind."* Python only. No C++ change, no new binding, no new dependency.

## What is missing

The engine produces an `IndexedMesh2` and nothing can save it. Every check on
the engine's output today is a triangle count, a status string or an SVG
picture of a hand-typed fixture. Nobody outside this repository can open
anything it makes.

```sh
git grep -nE "def write_|\.write_bytes|\.write_text" src_python/tin_engine
```

returns the SVG write in `cli.draw` and nothing else.

## Ruling 1 — the format is PLY

**Binary little-endian PLY, with an ASCII mode behind a flag.** The reasoning,
in the order it decided the question.

**QGIS reads meshes through MDAL, and MDAL does not read VTK or Gmsh.** VTK was
the first choice and it is wrong. MDAL's driver list is: 2DM, XMS TIN, NetCDF,
GRIB, XMDF, XDMF, DAT, 3Di, UGRID, FLO-2D, Selafin, HEC-RAS, SWW, Esri TIN,
SAGA FLOW, ADCIRC, **PLY**, DFSU, DFS2, H2i, MIKE 21. Neither `.vtu` nor `.msh`
is on it. Writing VTK would have produced a file the project's most likely
reader cannot open, and nobody would have found out until someone tried.

**PLY carries the edges.** MDAL's PLY driver reads and writes "any consistent
combination of vertices, faces and edges", with scalar and vector datasets on
each. That is the one property this engine actually needs from a format: the
constrained edges and their feature bits are the thing the whole pipeline
exists to preserve, and in most mesh formats they are the thing that gets
dropped. A PLY edge element with a scalar property per edge carries them.

**ParaView reads PLY too.** One format, both tools, no second writer.

**No new dependency.** PLY is a short ASCII header followed by packed records.
It is written with `numpy.ndarray.tobytes` and an f-string, in well under a
hundred lines. That matters here more than elsewhere: GDAL is prohibited by
`CLAUDE.md` §2, `meshio` would drag a format zoo in for one writer, and the
legacy route (below) needed `h5py` on top of that.

**Binary is the default** because ASCII costs both size and parse time, and the
parse time is the part that hurts. For a million-vertex TIN — so about two
million triangles, by Euler — the arithmetic is:

| | vertices | faces | total |
|---|---|---|---|
| binary, `double` x/y/z | 1e6 × 24 B = 24 MB | 2e6 × 13 B = 26 MB | ~50 MB |
| ASCII | ~40 B/line = 40 MB | ~25 B/line = 50 MB | ~90 MB |

Arithmetic, not a measurement. The user has ruled that no benchmark is a
precondition for this increment, and none is required by anything below: the
default is the format that reads as a block rather than as three million lines
of text, and the flag exists for the case where a person wants to read the file
by eye. A small mesh under `--ascii` is inspectable with `head`; that is the
whole reason the flag exists.

## Ruling 2 — coordinates are `double`, and float32 is not an option

**`property double x|y|z`.** float32 has about seven significant digits. At a
UTM 33N easting of 430 000 — `viz.fixtures.ORIGIN`'s, and representative of the
real data — 430 000 lies between 2¹⁸ and 2¹⁹, so consecutive float32 values
there are 2¹⁹ × 2⁻²⁴ = **3.1 cm** apart.

`cli.DEFAULT_SNAP_SPACING` is 1e-3 metres. So float32 is thirty times coarser
than the lattice the noder resolved the input on: two nodes the noder
deliberately kept apart can land on the same float32 coordinate, and a mesh
whose whole point is that crossings were resolved exactly would be saved with
the resolution thrown away. That is not a trade, it is a corruption.

The alternative, if a reader turns out to refuse `double`, is **not** float32
in world coordinates — it is a per-file origin offset: subtract a constant
easting and northing, write the residuals, and record the origin in a header
comment. `viz.fixtures.ORIGIN` is the precedent for the same idea in the other
direction. Residuals of a few kilometres in float32 hold millimetres
comfortably. This is written down as the fallback, not adopted: it costs a
header convention no reader understands automatically, and `double` is standard
PLY.

**Which makes one thing worth verifying before it is claimed.** Nothing in this
record has been opened in QGIS, because no session here can run QGIS. So: no
test and no document may say "QGIS opens this" until someone has opened one
(`PRINCIPLES.md` A1, A4). What the suite can assert is what the bytes are. What
a person must do once, by hand, is open a written file in QGIS and in ParaView
and say so here, with the version. If `double` is refused, the origin-offset
fallback above is the answer and it gets written up here with the failure that
forced it.

Byte order is forced little-endian by using `<`-prefixed numpy dtypes, so the
output does not depend on the host. A big-endian writer would be a second code
path with no reader here to test it against.

## Ruling 3 — two files, never one file holding both faces and edges

This is the caveat MDAL states about its own PLY driver: *"most host
applications (like QGIS) will expect the dataset to be either a 1d mesh — with
edges — or a 2D mesh — with faces."*

So a single PLY holding both is readable by MDAL and may still display as
nothing in QGIS. That is the worst available failure: valid file, silent
half-load, no error message.

**The ruling: each file holds one element type.**

- The surface file holds the vertex block and the face block. This is the 2D
  mesh, and it is what `--out` writes.
- The constraint file holds **the identical vertex block, in the identical
  order**, and an edge block instead of faces. This is the 1D mesh. It is
  written only when `--out-edges PATH` names a destination.

Because the vertex block is byte-identical in both, vertex *i* means the same
point in both files, and the two layers register on each other exactly when
loaded side by side in QGIS. That is a property of the design, not a
coincidence, and it is why the constraint file repeats the vertices rather than
referencing the other file — PLY has no cross-file reference, and inventing one
would be a private format.

Rejected alternatives:

- **One file with both, unconditionally.** Rejected on the caveat above.
- **One file, with a `--with-edges` flag putting edges in it.** Same defect,
  now opt-in. A flag whose effect is "this may silently not display" is worse
  than no flag.
- **Two files always.** Rejected: most callers want the surface, and a command
  that writes a file nobody asked for is a surprise. The second path is named
  or it is not written.

## Ruling 4 — z is data the caller supplies; this increment does not sample it

`IndexedMesh2` is 2D by design. Real elevations come from sampling the raster
per vertex, which is **MVP gap 3, "elevation assembly"**, and it does not
exist. Three ways out were on the table and one is right.

**The writer takes vertices as an `(N, 3)` float64 array and never invents a
z.** It does not take a callback, and it does not default to zero.

- Not a callback: that puts sampling policy inside a writer and makes the
  writer untestable without a raster. The arrow must point one way — the writer
  consumes an array, it never fetches one.
- Not a z=0 default: a default flat mesh is a wrong answer that looks like a
  right one, and it would be produced by the code path that has no idea it is
  wrong.

So increment 10 **does not depend on elevation assembly** and does not block on
it. When elevation assembly lands it supplies the array, and not one line of
the writer changes. That is the whole reason to draw the boundary at an array.

What this means for the CLI *today*, stated rather than discovered later: the
only meshes available to a command right now are the eleven gallery fixtures,
which are 2D shapes with no elevation source at all. So the command takes
`--flat`, required in the absence of any elevation, which fills z with 0.0 and
writes `comment elevation none (z=0, --flat)` into the header. Flat is
therefore always a word somebody typed and always visible in the file. The day
an elevation source exists, `--flat` becomes one of two ways to answer the
question and the command refuses to run with neither.

## Ruling 5 — CRS goes in a comment, and the file says it is not machine-read

PLY has no CRS field. The legacy wrote the projection as an HDF5 attribute
(`legacy/rasputin/tin_repository.py:169`,
`h5_points.attrs["projection"] = geometry.crs.to_proj4()`), which at least kept
it with the data.

Here the CRS is written as `comment crs <text>` in the header, and **no reader
acts on it**. QGIS will ask the user for the layer's CRS. That is a real
limitation of the format and it is recorded here so nobody spends a round
looking for the field. The comment is for the human who opens the file in six
months and needs to know what the numbers are in.

Where the text comes from: whatever the composition root holds, an EPSG code
preferred. The writer does not import `pyproj` and does not validate the
string — CRS handling is `raster.py`'s and `io/`'s job per
`project_structure.md`, and duplicating a fragment of it inside a byte writer
would be a second authority.

## Ruling 6 — where it lives, and what it is allowed to know

```
src_python/tin_engine/io/ply.py     # new. bytes in, bytes out, no path
src_python/tin_engine/cli.py        # the command, the only thing with a Path
```

`project_structure.md` marks `io/` as "all file decoding lives here (planned)".
This is the first encoding module in it, so that line is widened to decoding
**and encoding** in this PR. It stays the right home: `io/` is the layer that
knows file formats and knows nothing about `_core`.

The module is pure and imports numpy and the standard library only:

```
def write_ply(
    vertices,          # (N, 3) float64
    *,
    faces=None,        # (T, 3) uint32, or None
    edges=None,        # (E, 2) uint32, or None
    edge_properties=None,   # (E,) uint32, or None
    ascii=False,
    comments=(),
) -> bytes
```

Exactly one of `faces` and `edges` may be given — ruling 3, enforced at the one
place it can be enforced, with a `ValueError` naming both. `edge_properties`
without `edges` is the same error.

It returns `bytes` and takes no path, which keeps `06-cdt-viewer.md`'s "no file
is written below `cli.py`" true for the mesh path as it is for the SVG path,
and keeps the whole format testable with no filesystem, no raster and no
compiled extension: hand it three arrays and compare the bytes.

Data flow:

```
mesh + z array        -> np.column_stack        -> (N, 3) vertices
mesh.triangles                                  -> (T, 3) faces
mesh + noded PSLG     -> the constraint join    -> (E, 2) edges + (E,) masks
arrays                -> io.ply.write_ply       -> bytes      (pure)
bytes + Path          -> cli writes             (the only I/O)
```

The constraint join — which mesh edges are constrained, and which feature bits
each carries — already exists as `viz.scene.build_scene`, whose `SceneEdge`
carries `constrained`, `role` and `properties`. **It is not called from here.**
`viz/` is the renderer; a mesh writer reaching into it for a topology join
would couple output to visualization and make `viz/` load-bearing for a path
that draws nothing. The honest options are to lift the join into a shared pure
helper, or to compute the mesh-edge set here from `constrained_edges` alone and
carry only the mask. **This increment takes the second**, and the reason is
scope: the edge file's property scalar is worth having, and refactoring
`viz.scene` to share a join is a second increment's change to a module with a
live suite. `@developer` does not touch `viz/`.

Consequence, stated: the edges written are the mesh's own constrained edges,
per `IndexedMesh2.constrained_edges`, deduplicated across the two triangles
that share each one. Their feature bits come from the noded PSLG's
`edge_properties` looked up by node pair. If a constrained mesh edge has no
entry there, its property scalar is 0, which `_core.pyi` already defines as
*unclassified* rather than *wrong*.

## Prior art in `legacy/`

```sh
grep -rlniE "\.ply|\.off|\.obj|\.stl|vtk|gmsh|\.msh|write_mesh|export|save" legacy/
```

returns:

```
legacy/rasputin/tin_repository.py
legacy/rasputin/geo_tiff_reader.py
legacy/rasputin/application.py
legacy/rasputin/writer.py
legacy/tests/test_read_raster_file.py
legacy/tests/test_tin_repository.py
```

Two of them are the mesh writers. **Nothing is carried across, and three
things are carried across as negatives.**

1. `legacy/rasputin/writer.py` writes XDMF through `meshio`
   (`from meshio import XdmfTimeSeriesWriter, write_points_cells`). MDAL does
   read XDMF, so the format choice was not wrong — the dependency was. XDMF is
   XML plus HDF5, which is two files and two libraries where PLY is one file
   and none.
2. `legacy/rasputin/tin_repository.py` stores TINs as HDF5 through `h5py`,
   with a hand-built XDMF sidecar assembled by string surgery over
   `xml.etree`. It is a repository format, not an interchange format, and this
   increment is interchange.
3. The legacy passed the mesh as `triangulate_dem.point3_vector` — a bound C++
   type — straight into the writer. Here the writer takes numpy arrays and
   names no engine type, so it is testable with no extension in the process.

The shape worth taking, and the only one, is the *idea* that per-face and
per-vertex data ride along with the geometry. Legacy did it with `shades` and
land-cover fields; PLY does it with extra properties on the element. That is
what makes ruling 3's property scalar cheap.

## Files and LOC

**Estimate, not a measurement** (`PRINCIPLES.md` B4).

| File | What | Est. non-comment lines |
|---|---|---|
| `src_python/tin_engine/io/ply.py` | header emission, binary and ASCII bodies, the one-element-type check | ~75 |
| `src_python/tin_engine/io/__init__.py` | new package, one re-export | ~3 |
| `src_python/tin_engine/cli.py` | the `mesh` command: fixture, flags, two destinations | ~40 |
| `src_python/tin_engine/cli.py` | the constrained-edge set and its property lookup | ~20 |
| | **total** | **~138** |

About 20% of `CLAUDE.md` §2's ceiling. Measure it with the instrument the
other records use, noting its known blind spot on a bare `*,`:

```sh
git diff master...HEAD -- src_python/tin_engine \
  | grep '^+' | grep -v '^+++' | sed 's/^+//' \
  | grep -vcE '^\s*(//|#|\*|/\*|\*/|$)'
```

Tests are excluded. Estimated at ~180 lines, the largest single piece being the
binary round trip.

## Testing

**The invariant-critical suite is the binary body**, and it is the one place a
mutation round is worth paying for, per `docs/increments/README.md`'s cost
constraint. A wrong offset, a wrong dtype or a swapped index produces a file
that is structurally valid and geometrically wrong, which is exactly the defect
class a count-based assertion misses.

What the suite pins:

- **A round trip through a reader written in the test**, not through the
  writer's own logic. The test parses the header it was handed, unpacks the
  body by the header's own declarations, and compares arrays with
  `assert_array_equal`. A test that re-derives offsets the same way the writer
  did would agree with the writer's bug (`PRINCIPLES.md` B1).
- **ASCII and binary produce the same geometry.** Parse both, compare arrays.
  This is what catches one mode drifting when the other is edited.
- **`double` really is 8 bytes on the wire**, checked by writing a coordinate
  that float32 cannot represent — 430000.001 — and reading it back equal. That
  probe fails if anyone quietly switches to float32, which is ruling 2's whole
  content and is otherwise invisible.
- **Both element types together is refused**, with a message naming both.
- **The edge file's vertex block is byte-identical to the surface file's.** The
  two files' registration in ruling 3 rests on it, so it is asserted rather
  than assumed.
- **The header carries the CRS comment and the `--flat` comment** when those
  were supplied.

`--out` and `--out-edges` both go through the existing `cli._destination`
boundary, so its refusals are already covered and are not re-tested here.

## Acceptance

```sh
rasputin mesh catchment --flat --out /tmp/catchment.ply \
                                --out-edges /tmp/catchment-edges.ply
```

writes two files; a person opens both in QGIS as mesh layers and sees the
surface with the constraint lines over it. The engine's output leaves the
repository for the first time.

The QGIS half of that is the manual verification ruling 2 names. It is written
into acceptance rather than into a test because no CI runner here has QGIS.

## Not in scope

- **A PLY reader.** Nothing in this project consumes a mesh file.
- **Raster ingestion and elevation assembly.** MVP gaps 1 and 3, both still
  without records. Ruling 4 is what keeps this increment independent of them.
- **Per-triangle datasets** — slope, aspect, DEM error. PLY carries them and
  the format work here makes them cheap later, but nothing computes them yet.
- **Any second format.** VTK is settled against by ruling 1 and does not come
  back without a reader that needs it.
- **CRS validation or reprojection.** Ruling 5.
- **Any change to `_core`, the bindings, the stub file, or `viz/`.**
- **A refinement knob.** MVP gap 2.

## The record this updates

`ROADMAP.md` gains a row 10, in this PR, per `docs/increments/README.md`'s
merge rule, and its MVP-gap 4 entry is the one this increment closes.

`project_structure.md`'s `io/` line changes from decoding to decoding and
encoding, and gains `ply.py`.
