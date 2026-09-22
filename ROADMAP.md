# Rasputin roadmap

An index, and only an index. What each increment *is* lives in its own record
under `docs/increments/`, and nothing is restated here: a duplicated detail is
how the previous version of this file rotted into a seven-year-stale description
of the CGAL project this one replaced.

**Authoring order is not ship order.** Increment 6 is designed sixth and ships
before `05b`, because the viewer depends on nothing the noder produces and 5b is
the increment that most needs a picture to be checked against. The reasoning is
in `docs/increments/06-cdt-viewer.md`, "The ordering ruling".

| # | What it is | Status | Record |
|---|---|---|---|
| 1 | Exact geometric predicates | shipped (`d0946bc`, `aeb667a`) | `docs/increments/01-predicates.md` |
| 2 | 2D geometry value types | shipped (`1a60ad4`) | `docs/increments/02-core-geometry.md` |
| 3 | The PSLG and its validator | shipped (`d984ff8`) | `docs/increments/03-pslg.md` |
| 4 | The constrained Delaunay triangulation | shipped (`f22ddd3`) | `docs/increments/04-cdt.md` |
| 5a | The noder's numeric floor: snap grid and hot pixels | shipped (`e090909`, #69) | `docs/increments/05-noder.md` |
| 5b | The noder proper: segment splitting and the noded PSLG | shipped (`93fc772`, #76) | `docs/increments/05b-noder-driver.md` |
| 5c | Wiring the noder through: CDT signature, bindings, CLI | **in review** (#77) | `docs/increments/05c-noder-wiring.md` |
| 5d | The corner graze: a single-point cell touch stops counting | designed, **unscheduled** | `docs/increments/05d-corner-graze.md` |
| 6a | The pybind11 CDT surface | shipped (`311459b`, #71) | `docs/increments/06-cdt-viewer.md` |
| 6b-i | `build_scene`: the mesh-and-PSLG-to-geometry mapping | shipped (`ab4681c`, #72) | `docs/increments/06-cdt-viewer.md` |
| 6b-ii | The SVG renderer, the fixture gallery and `rasputin draw` | shipped (`4834568`, #74) | `docs/increments/06-cdt-viewer.md` |
| 7 | Edge property sets, replacing the one-bit `is_river` | shipped (`2e7577c`, #75) | `docs/increments/07-edge-properties.md` |
| — | `raster/`: grid-to-world geometry and bilinear sampling | shipped (`7785fea`), **no record** | none — predates the protocol |

## What stands between here and an operational MVP

Named, not described. None of it has an increment record, and that absence is the
status: `docs/increments/README.md` step 1 is where each of these starts.

An MVP is one command turning a DEM and a catchment polygon into a terrain TIN
file. Five things are missing, in dependency order.

1. **Raster ingestion, Python side.** The C++ `raster/` module samples; nothing
   decodes a GeoTIFF into it. `project_structure.md` names `raster.py` as the
   only adapter from decoded data into `_core`, and marks it planned. This is
   where CRS stops.
2. **Refinement.** The largest piece and the actual product: coarsen a dense DEM
   under an error budget instead of triangulating what you are given. Plan in
   `parallel_refinement.md`, which is older than every increment record here and
   has none of its own. The legacy `-ratio 0.4` was this knob.
3. **Elevation assembly.** `IndexedMesh2` is 2D by design and z comes from
   sampling the raster per vertex. Both halves exist; nothing joins them.
4. **Mesh output.** There is no writer of any kind.
5. **A CLI that does the job.** `rasputin` has `version` and `draw`; `draw`
   renders built-in fixtures only. No path from a file on disk to a mesh on disk.

Open and not MVP-blocking: clipping to the catchment polygon
(`auto_catchments.md`, and `legacy/rasputin/geometry.py`'s Shapely intersection
is the prior art), 5d above, and land-cover partitioning, whose foundation is
increment 7's property sets and whose consumer does not exist.

## Why this file was wrong, and what now keeps it honest

Until this edit the table said 5b was "not started; record not yet written" and
6b-ii was "in flight". Both had merged. 5c, 5d, 7 and the whole `raster/` module
were absent. Four increments shipped without this file noticing.

The cause is structural rather than careless: **nothing in the protocol touched
it.** `grep -n ROADMAP docs/increments/README.md .claude/REQUIRED-READING.md
CLAUDE.md` returned nothing before this commit. Every increment record is
enforced by the loop — design, red, green, review — and the index that points at
those records was enforced by nobody, so it decayed exactly as fast as the work
moved. `docs/increments/README.md`'s merge step now names it.

`docs/increments/README.md` is the protocol those records follow, and it — not
this table — is the file that states how work gets done.
