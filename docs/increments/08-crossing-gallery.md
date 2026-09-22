# Increment 8 — the crossing gallery

Status: design settled. No production code, no tests, no engine change and no
C++ line. This increment adds three fixtures to the gallery, renames one, adds
one stylesheet rule and its precedence entry, corrects two false sentences in
`06-cdt-viewer.md` and one false comment in `fixtures.py`.

Depends on 5b, 5c and 7, all shipped. Depends on nothing unshipped.

## Why this exists

The user's requirement, over two conversations: *"the case where a road crosses
a river, or enters into a forest"*, then *"roads entering forests, bridges over
lakes, structures entering or leaving a catchment"*.

5b and 5c made all of it work. The gallery does not show it. That is a product
gap, not a code gap. `06-cdt-viewer.md` states the gallery's purpose; this
increment is that purpose applied to the cases the engine was actually built
for.

The gap is measurable rather than asserted. Every fixture in the gallery today
is judged by running the engine over it and reading the scene back:

```sh
.venv/bin/python -c "
from tin_engine.viz.fixtures import GALLERY
from tin_engine import cli
from tin_engine.viz.scene import build_scene
for n, f in GALLERY.items():
    a = cli._triangulated(f, True, cli.DEFAULT_SNAP_SPACING)
    s = build_scene(a.source, a.mesh, ok=a.ok, closed_roles=a.closed_roles)
    print(n, s.kind.name, [x.kind.name for x in s.findings])
"
```

Read the output for two things this increment changes: whether any fixture
carries a chain whose role is a **closed breakline** (an area feature that is
not a hole), and whether any fixture produces a `Finding`. At the branch point
the answer to both is no.

## Ruling 1 — the terrain under a bridge is the water surface

**A bridge is not grade-separated in this engine, and must not be.**

A TIN is a height field over the plane. `z` comes from sampling the raster, and
over a lake the raster gives water height. A bridge polyline is a constraint on
the **terrain's** triangulation, not a structure carrying its own elevation, so
it nodes into the shoreline exactly like any other crossing. The deck height is
infrastructure. The mesh was never going to hold it, and a mesh that tried
would stop being a height field.

This is written down because the obvious reflex is the wrong one. The reflex is
a grade-separation mechanism — a flag or a rule marking some crossings as *must
not node*. That mechanism would need a third state in the noder, a per-chain
elevation source, and a definition of what the triangulation means when two
constraints occupy one planar point. None of it is needed, none of it is
designed, and none of it should be built for this reason.

The consequence for the design is that **`bridge-over-lake` requires no engine
change whatsoever**. It is `road-crosses-river` with the river closed into a
ring. That is the fixture's entire justification: it is the row that looks like
it should be special and is not, and the picture is the argument. A reader who
sees a bridge triangulate into the shoreline with triangles under the lake does
not go looking for a grade-separation feature.

If the deck height is ever wanted, it is a separate product — a structure layer
draped over the terrain — and it does not enter the CDT. That sentence is the
boundary, not a plan.

## Ruling 2 — an area feature is a closed breakline, not a hole

The gallery has three roles: `outer`, `hole`, `breakline`. A forest is not a
hole: the terrain inside it is meshed, and excluding it would leave a void in
the height field. Nor is a lake, under ruling 1.

So an **area feature is a closed breakline**: a chain of role `Breakline` whose
last index repeats its first. Three facts make this legal rather than a guess,
and each is checkable:

- `pslg_builder.hpp`'s stored-closure stage is **exempt for breaklines** — read
  the comment at the `StoredClosure` check and the `is_closed(ch.role)` guards
  around it (`grep -n 'StoredClosure\|is_closed' include/terrain/core/pslg_builder.hpp`).
  An outer or hole ring may not store its closing index; a breakline may.
- Winding is not checked for breaklines, by the same guard. A closed breakline
  may wind either way.
- The closure is written by **repeating the index**, not by repeating the
  coordinate. Reusing the index is exact. Repeating the coordinate relies on the
  noder's snap grid to merge two vertices that happen to coincide, which is a
  correctness argument where none is needed.

**The closure must be explicit or the ring is not a ring.** A four-index chain
`[4,5,6,7]` is an open polyline missing one edge, and the mesh differs: run the
probe under "Evidence" below with and without the repeated index and compare the
node counts. The missing edge is the one the crossing constraint would have cut.

`fixtures.py` therefore gains its first closed-breakline chains. `viz/scene.py`
needs nothing: `closed_roles` is `("outer", "hole")` for a fixture and
`(ChainRole.Outer, ChainRole.Hole)` for a core graph, and a closed breakline is
correctly **not** in either list — its closing edge is stored, so the scene must
not synthesise one. Confirm by reading `cli.FIXTURE_CLOSED_ROLES` and
`cli.CORE_CLOSED_ROLES`; the check that would catch a mistake here is
`test_cli_draw.py::test_a_well_formed_fixture_draws_no_disagreement`, which
would report every ring closure as a `MASKED_EDGE_WITHOUT_CHAIN`.

## Ruling 3 — `not-noded` is renamed to `road-crosses-river`

`06-cdt-viewer.md` ruled the name should stay, on the grounds that it describes
the input. That ruling predates 5c, which changed the fact it was made against.

The name is wrong for a simpler reason than being out of date: **it encodes a
resolved state of the tree.** Whether that input is noded is a property of which
commit you are standing on — it was false before 5b and is true after. `PRINCIPLES.md`
B2 is about citations, but the argument is the same one: a name that resolves to
a value has an invisible expiry date, and this one has already expired once.

`road-crosses-river` names the **geometry**, which no commit can change. It is
true whichever state the engine is in, and it is the user's own sentence.

A person scanning `rasputin draw --list` for a correct crossing must find it.
Today the line reads `not-noded: a road crossing a river: two constraints
meeting at a point neither names, noded` — the name says failure, the
description says it was fixed, and the reader skips the line.

### What the rename does not touch

`05b-noder-driver.md` and `05c-noder-wiring.md` name `not-noded` and are **not
edited**. They are dated records of what was true at their merge, and rewriting a
shipped record to agree with a later tree is the move `PRINCIPLES.md` B3 warns
about: the correction is itself a change and can carry the next defect. The
rename is findable from two places instead — this file, and the row in
`06-cdt-viewer.md`'s gallery table, which is the gallery's own record and is
being corrected in this PR anyway.

The set of files that mention the old name, for whoever does the edit:

```sh
grep -rln "not-noded" . --exclude-dir=legacy --exclude-dir=.git \
  --exclude-dir=build --exclude-dir=build-pyext --exclude-dir=.venv
```

Of those, `fixtures.py`, `tests/python/test_cli_draw.py` and
`tests/python/test_viz_svg.py` are renamed; the two C++ suites use the phrase for
the `NotNoded` status and are untouched; the three increment records are as
ruled above.

## The fixtures

Three new rows, all authored in the existing `_BOX_700_500` frame so the three
pictures are one controlled comparison rather than three unrelated shapes.
Coordinates are local metres; `ORIGIN` puts them at UTM 33N magnitudes, per the
fixture rule `06-cdt-viewer.md` states.

Property masks are **bare numbers**, because `viz/` imports no vocabulary. The
bits are `features.DEFAULT_VOCABULARY`'s; read them with
`grep -n 'EdgeProperty(name=' src_python/tin_engine/features.py` rather than
trusting this paragraph.

### `road-enters-forest`

A road crossing an area feature's boundary and stopping inside it.

| chain | role | mask | indices | points |
|---|---|---|---|---|
| 0 | `outer` | 0 | `range(4)` | `_BOX_700_500` |
| 1 | `breakline` | 0 | `[4,5,6,7,4]` | `(200,150) (500,150) (500,350) (200,350)` |
| 2 | `breakline` | road bit | `[8,9]` | `(60,250) (350,250)` |

The forest ring carries **mask 0**, not a forest bit, because the vocabulary has
no forest bit. Adding one is increment 7's territory, not this increment's, and
a fixture that invented a bit number would be asserting something the vocabulary
does not say. The ring draws in the plain breakline stroke.

Minimal for what it shows: one crossing, at `(200,250)`, and the road's far end
at `(350,250)` inside the ring. One crossing rather than two is the minimum for
"a breakline crosses an area boundary", and the interior endpoint is what makes
it *entering* rather than *passing through* — which is the difference between
this row and `bridge-over-lake`. The dangling interior endpoint is not new to
the gallery — `sliver-fan`'s rays already end inside the domain; read its
coordinates against its outer ring — so it is not claimed as a gain here.

### `wall-leaves-domain`

A structure leaving the catchment: a breakline crossing the **outer ring**.

| chain | role | mask | indices | points |
|---|---|---|---|---|
| 0 | `outer` | 0 | `range(4)` | `_BOX_700_500` |
| 1 | `breakline` | wall bit | `[4,5]` | `(350,250) (900,250)` |

One segment, one crossing of the outer ring at `(700,250)`, and one endpoint
outside the domain. That is the whole case. A four-point box rather than the
`catchment` outline, because the fixture rule is that a fixture must be
checkable by reading, and the crossing of an eleventh outline edge is not.

**This row is the first in the gallery to draw a `Finding`, and that is
intentional.** The exterior half of the wall is a chain edge the mesh does not
carry a constraint mask for, so `scene._findings` reports it as
`CHAIN_EDGE_WITHOUT_MASK` and `svg.py`'s `line.finding` draws it red and dashed.
For a segment outside the domain that is the **correct** answer, not a defect:
the engine is saying "your input claimed this edge and the mesh does not have
it", which is exactly what happens to a structure that leaves the catchment.

The finding overlay shipped in 6b-ii and no gallery picture has ever shown one.
`06-cdt-viewer.md`'s own standard — a presentation nobody has looked at is a
presentation that is wrong — applies to it. This row discharges that.

The row's description must say so, because a reader who sees red without an
explanation concludes the engine is broken. Proposed text, and the description
is the row's load-bearing part here: *"a structure leaving the catchment: a wall
crossing the outer ring. The exterior half is drawn as a finding, because the
mesh correctly does not contain it."*

### `bridge-over-lake`

Ruling 1, drawn.

| chain | role | mask | indices | points |
|---|---|---|---|---|
| 0 | `outer` | 0 | `range(4)` | `_BOX_700_500` |
| 1 | `breakline` | coastline bit | `[4,5,6,7,4]` | `(250,150) (450,150) (450,350) (250,350)` |
| 2 | `breakline` | road bit | `[8,9]` | `(120,250) (580,250)` |

Two crossings, at `(250,250)` and `(450,250)`, and triangles under the lake. The
road passes clean through, which is what separates this picture from
`road-enters-forest` — that one stops inside, this one spans.

The shoreline carries the **coastline** bit, which the vocabulary does have.
That is what makes this a lake rather than a generic ring, and it is why this
increment needs one stylesheet token (below).

`road` rather than a bridge-specific bit, because a bridge is a road. There is
no bridge bit and this increment does not invent one; ruling 1 is the reason the
engine does not need to know.

## The one stylesheet token

**`coastline` is added to `svg.py`'s `STYLESHEET` and to `cli.py`'s
`_PRECEDENCE`.** Everything else in this increment reuses existing tokens.

`_PRECEDENCE` becomes `("river", "coastline", "road")` — water over
infrastructure, which is the rule `cli.py` already states for why the list is
named by feature rather than derived from bit order.

The CSS rule: the same water hue as `line.river`, at a heavier stroke width. A
shoreline is the same substance as a river and should read as water; it is an
**area** boundary and should carry more weight than a line feature. No dash —
`line.finding` owns dashing, and a dashed water line on a correct picture would
read as the alarm it is not.

Check what exists before writing the rule rather than trusting this paragraph:

```sh
grep -n '^line\.' src_python/tin_engine/viz/svg.py
grep -n '_PRECEDENCE' src_python/tin_engine/cli.py
```

### Why this is worth one line in `cli.py`

The brief for this increment scoped the stylesheet but not `cli.py`. The two are
one change and cannot be separated: a token with a CSS rule and no precedence
entry is never emitted, and a precedence entry with no CSS rule draws
identically to the row above it. `06-cdt-viewer.md` records that second failure
mode as the reason `_PRECEDENCE` was held at one element, and states the growth
path — *"it grows when the stylesheet does"*. This is that growth, taken once.

The cost is honest and should be stated: the legend is derived from
`style.property_strokes` and is global, so **every** picture in the gallery gains
one legend row, not just the two that carry the bit. Three property rows is
still far short of the six that made `06-cdt-viewer.md` refuse to declare the
whole vocabulary.

The alternative — leaving the shoreline unstyled — was considered and rejected.
It makes `bridge-over-lake` and `road-enters-forest` the same picture in the same
colours, and it removes the thing the bridge row exists to show. A reader has to
be able to see that the ring is water before "the terrain under it is the water
surface" is an argument rather than a caption.

## The correction in passing

`fixtures.py`'s `CATCHMENT` comments call the two holes "Lake one" and "Lake
two". They are not lakes. Under ruling 1 a lake is meshed terrain at water
height, so a hole cannot be one — the comment contradicts a ruling in its own
file. What they might be is a creek intake or anything else excluded from the
catchment; **this increment does not model that and does not guess.** The
comments become "The first interior hole" / "The second interior hole", which
says what the geometry is and asserts nothing about the landscape.

The fixture's own `description` string is untouched at the geometry level but
carries the same false noun (*"two lake holes"*) and is corrected with it.

Per `PRINCIPLES.md` C3 this lands in this PR or it is not recorded.

## What 5b and 5c made false in `06-cdt-viewer.md`

Found while writing this design; fixed in this PR, same rule.

1. **The `not-noded` row of the gallery table** calls it *"a **deliberate
   failure fixture**, rendering the non-`Ok` presentation"*. It renders a mesh.
   Re-run the probe at the top of this file and read that fixture's line.
2. **"Three of the eight are failure presentations, and they are rows 4, 7 and
   8"**, and the sentence after it calling `not-noded` a backend refusal of
   `NotNoded`. Two of the eight are, and they are rows 4 and 8. The same
   paragraph's closing remark — *"three refusal pages out of eight is more than
   this gallery set out to have"* — was made true by 5c and by this increment,
   not by any action taken on that concern.
3. **The two-property fixture's stated blocker.** That section says
   `cli.py`'s `_PRECEDENCE` is `("river",)` *"one element, because `svg.py`'s
   `STYLESHEET` carries exactly one property rule"*. Both halves are false now;
   resolve them with the two greps under "The one stylesheet token". The fixture
   it describes — one chain carrying two bits at once, so first-match-wins is
   visible in a picture — **is still owed and is still not in the gallery**, and
   this increment deliberately does not add it. It is a different question
   (draw precedence) from this increment's (crossings), and folding it in would
   make the row set answer two things. The blocker sentence is corrected; the
   debt is restated, not discharged.
4. **The `catchment` row** says *"two lake holes and a braided breakline"*.
   Neither noun survives: see the correction above for "lake", and the chain is
   a single open three-point polyline, not a braid. Read `CATCHMENT`'s chain
   list in `fixtures.py`.

`fixtures.py`'s module docstring says *"eight shapes"*, *"Two of the eight"* and
*"``not-noded`` was a third until increment 5c and is now the showcase"*. All
three sentences move with the gallery and are rewritten here. The last one is
worth keeping in substance — the regression value of a fixture whose picture
changed between two commits is real — restated under the new name.

## Scope

- `src_python/tin_engine/viz/fixtures.py` — three fixtures, one rename, the
  `GALLERY` mapping, the module docstring, the catchment comments.
- `src_python/tin_engine/viz/svg.py` — one line in `STYLESHEET`.
- `src_python/tin_engine/cli.py` — one entry in `_PRECEDENCE`, and the comment
  above it that explains why it grew.
- `docs/increments/06-cdt-viewer.md` — the four corrections above.
- `docs/increments/08-crossing-gallery.md` — this file.
- `ROADMAP.md` — the row for increment 8, per `PRINCIPLES.md` C4.
- `tests/python/test_viz_svg.py`, `tests/python/test_cli_draw.py`.

Explicitly out of scope: any C++ file, `bindings/`, `_core.pyi`,
`src_python/tin_engine/features.py` (no new vocabulary bit), `viz/scene.py`,
`viz/style.py`, `viz/protocols.py`.

## What is worth testing

`@tester` owns this; the design's claims, in the form an assertion can hold.

**Mechanical, from the rename and the three new rows.** `GALLERY_NAMES` in both
suites goes from eight to eleven and loses `not-noded`.
`test_it_holds_exactly_the_eight_named_fixtures` is renamed with its count.
`CLASSIFIED_FIXTURES` in `test_viz_svg.py` gains three entries with their
reasons, and `test_no_fixture_carries_a_property_by_accident` is held against
the new map.

**The closed breakline is a ring.** For each of the two area fixtures, the
chain's first and last index are equal and the chain has role `breakline`. This
is the assertion that fails if someone "tidies away" the repeated index, and the
tidying is silent — the fixture still builds, still triangulates `Ok`, and draws
a picture missing one edge.

**Each crossing produces a constructed node.** For each new fixture, the noded
PSLG contains a vertex at the predicted intersection that the fixture's own
vertex array does not. The predicate must be the *position*, not a count: a
count assertion passes on a node constructed in the wrong place. The predicted
positions are `(200,250)`; `(700,250)`; `(250,250)` and `(450,250)`, plus
`ORIGIN`.

**The probe must be able to fail** (`PRINCIPLES.md` A3). The controlled negative
for the crossing assertion already exists in the suite:
`test_viz_svg.py::TestGalleryPredicates` self-tests `properly_cross` against a
shared endpoint. Extend that habit — assert that `road-enters-forest`'s road
endpoint at `(350,250)` is *not* a constructed node, so an oracle that called
every vertex constructed is red.

**`bridge-over-lake` meshes the lake.** At least one triangle whose **centroid**
is strictly inside the shoreline ring. This is ruling 1 as an assertion, and it
is the one test that fails if someone later adds grade separation. A triangle
count alone does not hold it — the mesh has triangles either way.

The centroid rather than the vertices, and this is not a stylistic choice: the
lake region has **no vertex strictly inside it**. Every vertex bounding those
triangles lies on the shoreline or on the two constructed crossing nodes, which
are also on it. A "all three vertices strictly inside" oracle is red on correct
output. Check before writing the assertion:

```sh
.venv/bin/python -c "
import numpy as np
from tin_engine.viz.fixtures import GALLERY, ORIGIN
from tin_engine import cli
a = cli._triangulated(GALLERY['bridge-over-lake'], True, cli.DEFAULT_SNAP_SPACING)
v = np.asarray(a.mesh.vertices) - ORIGIN
inner = (v[:,0] > 250) & (v[:,0] < 450) & (v[:,1] > 150) & (v[:,1] < 350)
print('vertices strictly inside the lake:', int(inner.sum()))
"
```

**`wall-leaves-domain` draws exactly one finding, and it is the exterior
segment.** Both halves: the count, and the endpoints. `scene.findings` gives
kind and vertex pair; resolve the pair through the noded PSLG's vertex array and
check the segment is the one from `(700,250)` to `(900,250)`. A count-only
assertion passes on a finding raised against the interior half, which would be a
real defect.

**The `coastline` stroke reaches the picture.** The pattern
`test_the_river_fixture_draws_its_property_stroke` already uses, applied to
`bridge-over-lake`, plus the negative that `road-enters-forest`'s ring — same
geometry class, mask 0 — carries no property token.

**Draw precedence is unchanged by the new entry.**
`test_cli_draw.py::test_water_is_drawn_over_infrastructure` pins river over
road; the new entry sits between them and must not reorder that pair.

Not an invariant-critical suite in `docs/increments/README.md`'s sense, so **no
mutation round.** This is declarative data and one CSS line; the failure mode is
a typo, and `README.md`'s cost constraint says to spend mutation where the
topology decisions are.

## Evidence

The probe every geometric claim above was taken with, before the claim was
written. It builds each proposed fixture from the coordinates in the tables,
runs it through the real engine and prints what came back. It is re-runnable and
it is how `@tester` should check the tables before writing against them.

```sh
.venv/bin/python - <<'PY'
import numpy as np
from tin_engine.viz.fixtures import _fixture, ORIGIN, _BOX_700_500 as BOX
from tin_engine import cli
from tin_engine.viz.scene import build_scene

CASES = {
 "road-enters-forest": (BOX + [[200.,150.],[500.,150.],[500.,350.],[200.,350.],
                               [60.,250.],[350.,250.]],
   [(range(4),"outer",0), ([4,5,6,7,4],"breakline",0), ([8,9],"breakline",2)]),
 "wall-leaves-domain": (BOX + [[350.,250.],[900.,250.]],
   [(range(4),"outer",0), ([4,5],"breakline",32)]),
 "bridge-over-lake": (BOX + [[250.,150.],[450.,150.],[450.,350.],[250.,350.],
                             [120.,250.],[580.,250.]],
   [(range(4),"outer",0), ([4,5,6,7,4],"breakline",8), ([8,9],"breakline",2)]),
}
for name, (pts, chains) in CASES.items():
    f = _fixture(name, "probe", pts, chains)
    a = cli._triangulated(f, True, cli.DEFAULT_SNAP_SPACING)
    s = build_scene(a.source, a.mesh, ok=a.ok, closed_roles=a.closed_roles)
    given = {tuple(np.round(p, 6)) for p in np.array(pts) + ORIGIN}
    built = [tuple(np.round(p, 6) - ORIGIN)
             for p in np.asarray(a.source.vertices) if tuple(np.round(p, 6)) not in given]
    print(name, a.status, s.kind.name, "constructed:", built,
          "findings:", [x.kind.name for x in s.findings])
PY
```

The masks in the probe are bare numbers because `viz/` holds no vocabulary; they
must agree with `features.DEFAULT_VOCABULARY`, and the grep under "The fixtures"
is what checks that. The probe is *not* a substitute for the suite: it shares
`cli._triangulated` with the code under test, so it establishes that the
geometry is buildable, not that the renderer is right.

Two things the probe establishes that the design leans on and that are easy to
assume the other way:

- a closed breakline written with a **repeated index** builds, triangulates `Ok`
  and is exempt from the stored-closure check;
- dropping the repeated index also triangulates `Ok`, with a different mesh.
  Both pass. Only the assertion under "What is worth testing" separates them.

## Prior art in `legacy/`

**Nothing is carried across.** The commands and what they returned:

```sh
grep -rln "breakline\|constraint" legacy/
# legacy/rasputin/land_cover_repository.py
# legacy/rasputin/globcov_repository.py
# legacy/rasputin/triangulate_dem.h
# legacy/rasputin/wfs_repository.py
# legacy/rasputin/gml_repository.py

grep -rln "forest\|bridge\|lake" legacy/
# legacy/bindings.cpp
# legacy/rasputin/globcov_repository.py
# legacy/rasputin/material_specification.py
# legacy/rasputin/triangulate_dem.h
# legacy/rasputin/web_visualize.py
# legacy/rasputin/gml_repository.py

grep -rn "bridge" legacy/
# (no output)
```

What that set contains, read rather than inferred:

- **Forest is a raster land-cover class, not a vector area.**
  `legacy/rasputin/globcov_repository.py` enumerates `forest_type_1` through
  `forest_type_6` as GlobCover codes, and `land_cover_repository.py` maps a code
  to a colour, a description and a rendering material. A face got an attribute
  from a raster lookup. Nothing ever became a constraint chain, so nothing in
  the legacy tree has an opinion about a road crossing a forest boundary.
- **Lake is a rendering material.** `material_specification.py` defines
  `lake_material`. That is a shading dictionary for the web viewer, not a
  geometric feature, and it has nothing to say about ruling 1.
- **No bridge anywhere.** The grep returns nothing. The grade-separation
  question this file rules on was never posed in the legacy tree, which is
  consistent with the ruling — it does not arise until constraints cross.

The legacy per-face attribute approach is orthogonal to this increment rather
than superseded by it: an area *boundary* as a constraint and an area *class* as
a face attribute are two different mechanisms and a finished product will want
both. Nothing here forecloses the second.

## LOC

**Estimate, not a measurement** (`PRINCIPLES.md` B4).

| File | What | Est. non-comment lines |
|---|---|---|
| `viz/fixtures.py` | three `_fixture` calls, three `GALLERY` entries, one rename | ~40 |
| `viz/fixtures.py` | module docstring rewrite | ~12 |
| `viz/svg.py` | one `STYLESHEET` rule | 1 |
| `cli.py` | one `_PRECEDENCE` entry | 1 |
| | **total** | **~55** |

Comfortably under `CLAUDE.md` §2's ceiling, with no seam to pre-declare. The
estimate assumes the three fixtures reuse `_BOX_700_500`; if they do not, add
about nine lines.

Measure it against the ceiling with the instrument every figure in
`06-cdt-viewer.md` was taken with — noting, as that file does, that the
instrument counts docstring lines and has a known blind spot on a bare `*,`:

```sh
git diff master...HEAD -- src_python/tin_engine/viz/fixtures.py \
    src_python/tin_engine/viz/svg.py src_python/tin_engine/cli.py \
  | grep '^+' | grep -v '^+++' | sed 's/^+//' \
  | grep -vcE '^\s*(//|#|\*|/\*|\*/|$)'
```

Tests are excluded from the ceiling. Estimated at ~120 lines across the two
suites, most of it the mechanical `GALLERY_NAMES` widening.

## Open, and deliberately not closed here

- **The two-property fixture** `07-edge-properties.md` asked for. Unblocked as
  of this increment's stylesheet change, still not in the gallery, and still a
  question about draw precedence rather than about crossings. See correction 3.
- **A `forest` bit in `DEFAULT_VOCABULARY`.** The forest ring draws at mask 0
  because no honest number exists. Whether the vocabulary should carry area
  classes at all is increment 7's question reopened, not this one.
- **`hole-in-hole` as a mesh.** `06-cdt-viewer.md` records that the design meant
  a mesh and the red step narrowed it to a refusal, and parks a ninth fixture or
  a widened test as a follow-up. This increment does not take it: it is a
  question about in-domain nesting, not about crossings, and taking it would
  make the row set answer two things.
