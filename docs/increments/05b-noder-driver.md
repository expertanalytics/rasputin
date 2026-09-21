# Increment 5b — the noder's topology half

Status: design settled and **amended after the red round**, which refuted one of
its own claims. `@tester`'s suite found guarantee 14 false on this document's
own worked degeneracy; the section that asked for the check is corrected in
place rather than reworded, the consequence is ruled on as risk 16, and two
further findings — `RingDegenerateAfterSnap`'s reachability and
`NodedPslgBuilder::build`'s status for a guarantee-14 refusal — are settled in
the sections that own them, the first of them by a bounded exhaustive search
whose program and output are inline. Three smaller measured corrections are folded into
mutants 11 and 12 and the builder's suite. No production line count moves.
Also: **5b is split before the red step, not after it.**
`docs/increments/05-noder.md` designed 5a in full and fixed 5b's seam; this
document discharges that obligation and overturns three things in it on
evidence. The re-estimate below comes to **~786 non-comment production lines**
against `CLAUDE.md` §2's 700, so the contingency split `05-noder.md` names is
taken here, at design time — the last moment at which cutting a seam is cheap,
for the reason "The LOC gate" gives. The split point is **not** the one
`05-noder.md` proposed, and the reason is measured rather than argued.

The two PRs:

- **5b — the noder, in C++.** `core/noded_pslg.hpp`, `noding/broad_phase.hpp`,
  `noding/noded_pslg_builder.hpp`, `noding/node.hpp`. **Purely additive**: no
  existing header changes, nothing in `src/`, no binding moves, the build stays
  green at every commit.
- **5c — the crossing.** `cdt::triangulate`'s signature changes from
  `const Pslg&` to `const NodedPslg&`, and the noder crosses into Python:
  `bindings/core.cpp`, `src_python/tin_engine/_core.pyi`,
  `src_python/tin_engine/cli.py`. This is the PR in which
  `rasputin draw not-noded` stops drawing zero triangles.

Three premises from `05-noder.md` are overturned below, each with the command
that refutes it:

1. **The signature change is not mechanical and is not ~10 lines.**
   `bindings/core.cpp:466` calls `terrain::cdt::triangulate<DetriaBackend>(pslg,
   options)` on a `const Pslg&`. Retyping the CDT entry point breaks that
   translation unit, and the only repair is to bind a producer of `NodedPslg` —
   which drags `_core.pyi`, `cli.py`, the Python suites and the renderer's scene
   join along with it. `05-noder.md`'s "Files and LOC" table books this as one
   row, `include/terrain/cdt/*` plus `src/cdt/detria_backend.cpp`, ~10, and
   names no Python file at all. See "What the signature change actually costs".
2. **Guarantee 14(a) as written is false on the one input guarantee 15 exists
   for.** It requires `classify<K>` to return `Disjoint` or `Touching` for every
   candidate pair; a road noded along a river yields two chains carrying the
   *same* edge, and `classify` on two identical segments returns `Overlapping`.
   The amended clause is below and it is strictly more useful: duplicate edges
   are permitted, partial overlaps are not, and after splitting there is no
   third case.
3. **The project has two LOC instruments and they disagree by 15 % on C++.**
   `docs/increments/06-cdt-viewer.md` calls
   `grep -vcE '^\s*(//|#|\*|/\*|\*/|$)'` "the instrument every figure in this
   document was taken with"; its `#` alternative is a Python comment and a C++
   *preprocessor directive*, so on C++ it silently drops every `#include` and
   `#pragma once`. Measured on 5a: 153 by that command against the 179 the green
   commit `41ac5b6` reports, and
   `for f in include/terrain/core/snap_grid.hpp include/terrain/noding/intersect.hpp include/terrain/noding/node_set.hpp; do grep -cE '^\s*#' $f; done`
   prints `6`, `10`, `10` — exactly the 26-line difference. Instruments are
   declared per language under "Files and LOC" and the figures below are not
   interchangeable with 6b's.

## Prior art in `legacy/`

**Nothing, and the evidence is that the legacy pipeline never noded at all.**

```sh
grep -rln "noding\|snap\|intersect\|broad phase\|split" legacy/
```

returns **six** files: `legacy/bindings.cpp`,
`legacy/rasputin/triangulate_dem.h`, `legacy/rasputin/reader.py`,
`legacy/rasputin/geometry.py`, `legacy/rasputin/gml_repository.py`,
`legacy/tests/test_polygons.py`. Every hit is one of three things, and none is a
noder:

- **Polygon boolean operations delegated to a library.** `CGAL::intersection` at
  `legacy/bindings.cpp:94` and `legacy/rasputin/triangulate_dem.h:414`; Shapely's
  `.intersection` / `.intersects` at `legacy/rasputin/geometry.py:243-251`,
  `legacy/rasputin/gml_repository.py:175-178` and
  `legacy/rasputin/reader.py:434-468`, and the test of it at
  `legacy/tests/test_polygons.py:25`. These clip a domain against a polygon.
  They compute no arrangement, split no constraint and produce no node.
- **Ray–mesh intersection for shadows**, `legacy/rasputin/triangulate_dem.h:614`.
  A different question entirely.
- **String and attribute splitting.** `.split()` on GML coordinate text at
  `legacy/rasputin/gml_repository.py:147`, and `split_by_colors` at
  `legacy/rasputin/geometry.py:46`, which partitions faces by a colour
  attribute.

`grep -rni "snap" legacy/` and
`grep -rni "noding\|broad.phase\|hot.pixel" legacy/` both return **nothing** —
zero lines, in a tree of 1,422 lines of domain logic.

The positive finding is more useful than the negative one.
`legacy/rasputin/triangulate_dem.h:49` declares

```cpp
using ConstrainedDelaunay = Constrained_Delaunay_triangulation_2<Gt, CGAL::Default, CGAL::Exact_predicates_tag>;
```

and `:485` feeds it `dtin.insert_constraint(...)`. `Exact_predicates_tag` is the
CGAL tag under which the CDT *itself* detects intersecting constraints and
inserts the intersection points. **Noding is precisely the capability that left
with CGAL**, and the legacy tree contains one line of it — a template argument.
There is nothing for `@migration-expert` to read and nothing to carry across;
`@migration-expert` is not spawned for 5b or 5c.

That is also why the failure the increment fixes is *new*. It was never a
failure in the legacy product, because the library absorbed it silently — which
is the same silence `04-cdt.md`'s ruling refuses ("never silently produces a
mesh that ignores a crossing").

## The target case, verified on this tree

Master at `d6c7c88`, the built extension in `.venv`:

```sh
.venv/bin/python -c "
import numpy as np
from tin_engine import _core as c
V=np.array([[0,0],[700,0],[700,700],[0,700],[100,100],[600,600],[100,600],[600,100]],dtype=float)+np.array([430000.,6900000.])
ch=[(list(range(4)),c.ChainRole.Outer,False),([4,5],c.ChainRole.Breakline,False),([6,7],c.ChainRole.Breakline,False)]
r=c.build_pslg(V,ch); print('build ok:', r.ok)
o=c.triangulate(r.pslg); print(o.status, o.mesh.triangle_count, o.message)
"
```

prints

```
build ok: True
CdtStatus.NotNoded 0 Constrained edges are intersecting: first edge between points 6 and 7, second edge between points 4 and 5
```

Those are the coordinates of the `not-noded` gallery fixture
(`src_python/tin_engine/viz/fixtures.py`, on branch
`increment6b-ii-renderer-red` at `4482649`): a square outline and two breaklines
crossing at local `(350, 350)`, a point neither chain names. The `Pslg` validator
accepts it — correctly, `03-pslg.md`'s ruling is that a valid `Pslg` promises no
pairwise disjointness — and the CDT refuses it.

**The acceptance criterion for 5b+5c is that the same command, after the
signature change, returns `Ok` with a mesh whose vertex array contains the
snapped image of `(350, 350) + ORIGIN` and contains it once.** 5b establishes it
in C++; 5c is where the command above can be run.

## What the renderer buys this increment

`rasputin draw <fixture> --out x.svg` exists on branch
`increment6b-ii-renderer-red` (`4482649`) and is **not yet on master** — 5c
depends on 6b-ii merging, 5b does not depend on it at all. Three things it buys,
in descending order of how load-bearing they are:

1. **It has already falsified a piece of this design, before any code was
   written.** `src_python/tin_engine/viz/scene.py:200`'s `build_scene(pslg,
   mesh, ...)` joins the PSLG's chain edges to the mesh's masked edges **by
   vertex-index pair** (`_chain_edges` builds `dict[Pair, ...]` from
   `indices_of(c)`; `_mesh_edges` builds `set[Pair]` from triangle corners).
   After noding, the mesh's indices are **node ids** — `NodeSet` sorts by
   `GridPoint`, so they are neither the input vertex ids nor a permutation with
   any relation to them. Feeding the *input* `Pslg` to `build_scene` alongside a
   *noded* mesh therefore produces a join over two unrelated index spaces: every
   chain edge becomes a `CHAIN_EDGE_WITHOUT_MASK` finding and every constrained
   mesh edge a `MASKED_EDGE_WITHOUT_CHAIN`, and the whole gallery renders in the
   alarm colour with a correct mesh underneath it. **The scene must be built
   from the `NodedPslg`.** That is a 5c requirement, it is the reason
   `NodedPslg` must satisfy `viz/protocols.py`'s `PslgLike` structurally, and
   nothing in `05-noder.md`, `03-pslg.md` or `06-cdt-viewer.md` records it.
   `NodedPslg` having the same shape as `Pslg` was ruled for a different reason
   (the CDT wrapper); this is a second, independent one, and it is the one with
   a test behind it — `tests/python/test_viz_protocols.py` is the join.
2. **An unsplit T-junction is a line ending in the middle of another line**, and
   a person sees that in the time it takes the SVG to open. It is the one defect
   class in this increment with no cheap oracle: guarantee 14(b) is checked
   brute-force on a few dozen segments, and the gallery is where a survivor on a
   realistic shape would show up. `05-noder.md`'s own argument for building the
   viewer before 5b was that 5b is the increment that most needs an instrument;
   this is that argument cashed.
3. **`not-noded` becomes a regression fixture with two states.** Today it draws
   the PSLG alone with the backend's words in the header band
   (`SceneKind.FAILED`). After 5c it draws a mesh. The picture changing is the
   increment's product, and it is checkable by a person in one command rather
   than by reading an assertion count.

What the renderer does **not** buy: any relaxation of the property suites. A
right-looking wrong picture is exactly the failure mode `06-cdt-viewer.md` says
testing effort goes to, and the noder's silent failures — a dropped broad-phase
candidate, a lost edge property — are invisible at gallery scale.

## Edge properties: the type is increment 7's, the per-edge array is 5b's

The carrier is `terrain::EdgeProperties`, shipped by increment 7:
`include/terrain/core/edge_properties.hpp` for the type,
`src_python/tin_engine/features.py` for the vocabulary that names its bits, and
`docs/increments/07-edge-properties.md` for every ruling about either — the
width and its ceiling, the absence of named members, union rather than priority,
what was rejected, and the boundary model. **None of that is restated here.**
`docs/increments/README.md` asks a prompt to say "read the file" rather than to
carry the file; the same holds between two designs, and a copy is the thing that
goes stale. `git log --oneline -1 -- include/terrain/core/edge_properties.hpp`
names the commit that last moved the type.

What 5b consumes, unchanged:

- **`Chain::properties`**, an `EdgeProperties` — increment 3's `Chain::is_river`
  no longer exists. The noder preserves it chain-for-chain (guarantee 16).
- **Union as the merge**, commutative, associative and idempotent. Those three
  are the only properties of the type this increment leans on.
- **A core that names no feature.** There is no `River` enumerator in
  `terrain::`, so the noder cannot condition on one — which is why every ruling
  below is about *where* bits are merged and never about *which*.

**The noder does not care how many bits there are.** That was the test of whether
this design had been over-fitted to one bit, back when `Chain::is_river` was
still the field; it had not. The width change cost the sections below nothing:
same array, same index alignment against the flat `(c, k)` enumeration, same
`edge_base`, same dedup key, same one-directional oracle. Guarantee 14(a) — the
`Overlapping` amendment — is likewise untouched by it, being a statement about
geometry rather than about meaning.

Four things are 5b's own, and are ruled here. A fifth, deferred by this
design to whoever owned the stylesheet, has since been answered elsewhere and
closes the section.

### The array is dense: one entry per output edge, not a sparse override set

`edge_properties()` is a `std::span<const EdgeProperties>` of length
`edge_base(chains().size())`, one entry per flat `(c, k)`, which is guarantee 15.
`docs/increments/03-pslg.md` and `project_structure.md` each assumed a **sparse
override set** before this design and each now says dense, naming this document
for the shape.

An override set is sparse only if most edges agree with a single source value
that the exceptions override. **After noding there is no such source value.** An
output edge can descend from two input chains at once — that is the fact the
property *set* exists for — so the thing a sparse form would record an exception
*to* does not exist. Encoding an array that has no default sparsely buys nothing
and costs three things: an index structure, a branch on every read, and a second
spelling of "empty" that cannot be told apart from "absent". The empty set is
already the default `EdgeProperties{}`, and it means *unclassified* — a legal,
expected state, not a missing one.

Dense also costs nothing worth naming: one 32-bit word per output edge, in one
contiguous `std::vector`, trivially parallel-reducible, against an output that
already carries a `Point2` per vertex and a `std::uint32_t` per chain position.

### The merge happens at the node-id edge-key dedup, not at a classification

Step 5 of the split pass is where the union is taken: for each flat position
`(c, k)`, the union over every input chain carrying an output edge with the same
node-id pair `(min(u, v), max(u, v))`. The ruling's content is what it excludes:

- **No collinear-overlap classification participates.** The tempting spelling is
  "union the properties of every pair `classify<K>` calls `Overlapping`", and it
  would put a kernel predicate, and therefore floating point, inside the merge.
  It is also redundant: guarantee 14(a) as amended says two output edges that
  overlap at all have **equal** node-id pairs, so on noded output the integer key
  *is* the overlap relation. A predicate would re-derive, less exactly and at a
  cost, what the key already decides.
- **The key is integers**, which is `05-noder.md`'s standing rule for the merge,
  and is what keeps the result independent of which way a snap rounded.
- **The reduce is over an unordered collection.** Contributing chains arrive in
  chain-enumeration order and, within a round, in whatever order the broad phase
  populated buckets; the answer may not depend on either. Union's three algebraic
  properties are the whole licence for that, and they are why a priority scheme
  would have been unusable here even had one been on offer.
- **The dedup does not merge chains.** It is over the edge *set* and its product
  is the property array. Two chains that run together after snapping remain two
  chains, index-for-index with the input (guarantee 16); what they come to share
  is node ids, and therefore per-edge entries, not identity. The Python
  renderer already performs the same join one layer up — `viz/scene.py`'s
  `_chain_edges` unions `properties` per undirected pair while letting the first
  contributing chain keep the role — which is corroboration that the shape is
  right, not a dependency.

### Guarantee 15's oracle: built from the input, asserted as a subset

A design constraint on the suite rather than a choice the suite makes, so it is
ruled here; the assertion is written out under "What is worth testing", and the
reason 15 is not the builder's is under `noding/noded_pslg_builder.hpp`.

- **Built from the input, never from the dedup's own provenance map.** An oracle
  that asks the noder which chains it recorded as contributing, and then asserts
  the union over those, is the dedup restating itself: it is green on any merge
  that is internally consistent, including one that silently drops a contributor.
  Borrow the producer's *predicate* — `segment_meets_cell<K>` plus arc-order
  betweenness — and never its *records* (`computational-geometry/SKILL.md`; the
  worked history is `05-noder.md`, guarantees 14 and 15).
- **A subset, not an equality, and one direction only.** The union over the input
  chains that the relation admits for an output edge `e` is a subset of
  `edge_properties()[e]`. The converse is *false*, not merely unproven: a chain
  that passes near both nodes and spans them satisfies the relation without
  having contributed geometry, so the oracle's set can be strictly smaller than
  the truth. Asserting equality would make the suite red on correct output —
  increment 5's exact-incidence failure in a new costume.
- **And that is not the only reason. The oracle is exact only in round 1.** It
  splits the **input** segments against the **final** node set. The driver splits
  the input segments in round 1, but from round 2 it splits the **snapped pieces**
  — which step 6 mandates, and which is the whole of what makes 14(b) reachable
  at a fixpoint instead of merely re-detected. A piece is a different segment from
  the segment it came from: it has different endpoints, so a different arc order,
  a different `segment_meets_cell` answer for a node near the seam, and therefore
  a legitimately different refinement. On any input that takes a second round the
  oracle and the driver are relating different objects, and only the subset
  direction survives that. The first reason is about which *chains* the relation
  admits; this one is about which *segments* each side is splitting. Both point
  the same way, which is why one assertion covers them.
- **The converse's state, as evidence rather than as proof.** `@developer` could
  not construct an input on which the oracle claims a contribution the driver does
  not make, and no generated seed produced one. That is the state of the search,
  not a theorem: nothing here rules out a multi-round input where a piece's
  refinement drops a pair the input segment's would have had. If such an input
  exists the subset direction still holds, which is the point of choosing it.
- **This is the assertion to re-check when 5d changes the predicate.** Both
  reasons above are stated in terms of `segment_meets_cell<K>`: the first through
  "passes near both nodes", the second through what a piece's endpoints do to that
  answer. 5d narrows the predicate, so it moves the boundary that both arguments
  sit on, and a narrower predicate is exactly the direction that could turn an
  inert difference into a `proven` superset. Re-run the reasoning, not just the
  suite — the suite is green under both.
- A spuriously *added* property therefore escapes this assertion by construction
  and is owed a mutant instead; that is mutant 7.

### No bare word per edge, and the reason is now what the type is

5b may not carry an edge's properties as a `std::uint8_t`, a `bool`, or a
`std::uint32_t` anywhere inside `terrain::` — not on `NodedPslg`, not in the
builder's arrays, not as a span. When this prohibition was first written the
argument was that a one-bit `edge_is_river` would have to be rewritten once the
set arrived. That argument has expired. The stronger one is the shipped type:

**`EdgeProperties` has no conversion to or from `bool`, `int` or `std::uint32_t`
in either direction.** `bits()` is the one way out, `bit()` and `operator|` the
one way in, and the constructor that takes a word is private and two-argument so
that not even `std::is_constructible_v` finds a route. The header says why, and
the two failure modes it names are exactly the ones a bare per-edge array
reopens: `if (edge_is_river)` compiles again, and a call site that passed `true`
keeps compiling with a changed meaning. A bare array would reinstate both one
level below the field that closed them, on a per-edge array rather than a
per-chain field, and at a place where nothing re-validates: the mask admission
check lives in `bindings/core.cpp`, where untrusted masks arrive, and noded
output never passes through it.

Bare masks are legal in exactly two places, both outside `terrain::` — the
binding's admission of an untrusted mask, and the accessors that hand a mask to
Python. Between those two boundaries an edge's properties travel as the type.

### The renderer question this design deferred has been answered

An edge that is both a road and a river needed either a precedence order or two
strokes, and this document declined to choose, naming the decision as the
stylesheet owner's. Increment 7 made it: one stroke per edge under a declared
precedence — the order of `SvgStyle.property_strokes`, expressly not bit order,
since a bit position is a vocabulary's numbering and priority is the
stylesheet's — with the full set still carried on the edge as
`SceneEdge.properties`. See `docs/increments/07-edge-properties.md`,
"The renderer: the ruling 5b named and deferred". Nothing of it is 5b's, and
what remains 5b-adjacent is risk 15: at 5c the scene join must be fed the
`NodedPslg` and not the `Pslg`.

## The seam, and why it is not `05-noder.md`'s

`05-noder.md`'s contingency split was **5b** = `core/noded_pslg.hpp` +
`noding/broad_phase.hpp` + a `NodedPslgBuilder`, **5c** = `noding/node.hpp` +
the driver + the CDT signature change. It is rejected, on three grounds:

1. **It puts the signature change in the same PR as the driver**, and the
   signature change is the one that detonates the Python surface (premise 1
   above). That PR is the driver *plus* `bindings/core.cpp` *plus* `_core.pyi`
   *plus* `cli.py`, which is the largest of the three pieces, not the smaller
   half of a split.
2. **The measured density asymmetry runs along the language boundary, not along
   the type/algorithm boundary.** C++ increments have come in at or under
   estimate — increment 3's `pslg_builder.hpp` 249 against ~255, increment 4 267
   against ~342, increment 5a 179 against ~270. Python and binding work has come
   in at **twice** estimate: 6a shipped 467 against ~235, with
   `bindings/core.cpp` +283 against ~150 and `_core.pyi` +125 against ~60
   (`docs/increments/06-cdt-viewer.md`, "6a as shipped"); `scene.py` shipped 194
   against ~90, factor 2.16. A seam should be cut where the estimate is least
   trustworthy, and that is the C++/Python line.
3. **It leaves `NodedPslg` with no producer**, which `05-noder.md` itself flags
   ("`NodedPslg` would then need a producer to exist at all"). The seam below
   does not: 5b ships a complete, exercised noder.

The cost of the seam chosen here is the one 5a already paid and it is stated
rather than hidden: **5b ships with no production caller.** Only tests
instantiate `node<K>` until 5c merges. That is the second increment in a row
with that shape, and it is the direct consequence of the ceiling being per PR.

**The seam has a mechanism, like 5a's:** no file in 5b modifies an existing
file. `git diff --stat master...HEAD` on 5b's branch must show additions only,
under `include/terrain/` and `tests/`, plus `tests/cpp/CMakeLists.txt` and this
document. If `bindings/core.cpp` or `src/cdt/detria_backend.cpp` appears in that
list, the seam has been violated and the diff says so.

## The LOC gate, and the moment it fires

`docs/increments/06-cdt-viewer.md`'s second gate fired at 709 against 500, on
branch `increment6b-ii-renderer-red`, and **the split it authorised was not
taken**. Its own conclusion is the right one — "a LOC gate has to sit before the
red step, not after it... Any future pre-declared seam should be measured against
the *design's* file list at the moment `@tester` is briefed" — but the reason it
gives for it is wrong, and the wrong reason matters more than the right
conclusion, because it is what a later reader would act on.

**The reason given is that the split was *impossible*: one committed test file
(`test_viz_svg.py`) covered `style.py`, `svg.py` and `fixtures.py`, and the green
step may not touch a test file. `docs/increments/README.md:57-59` says the
opposite in so many words:**

> The rule is not "tests are frozen after red"; it is "`@developer` does not edit
> tests, and no test change hides inside an implementation commit".

Amendments after red are `@tester`'s, and land as their own commit with the
reason in the message — increment 3 pinned three behaviours that way and
increment 2 retuned three constants after implementation. So splitting
`test_viz_svg.py` was **available**: a `@tester` commit splitting the suite, a
second red round, a second PR. It was expensive, not forbidden.

**The honest generalisation is a cost argument, and it is stronger than the
impossibility one because it predicts what actually happened:**

> A LOC gate must sit before the red step, because the red step fixes the PR's
> shape. A gate placed after it is not unenforceable — the remedy is a `@tester`
> amendment splitting the suite, a second red round and a second PR — but that
> remedy costs more than the split it buys, so the gate will be argued away
> rather than obeyed.

Which is exactly the observed record. Three seams have been armed in this
project: increment 3's was named and never used; 6b-i's fired and was obeyed;
6b-ii's fired and was not. **One for three, and the one that worked fired before
its suite existed.** That record is the argument for the shape of the gates
below, and it is also an argument against arming many of them.

The impossibility framing lives in a paragraph that exists only on
`increment6b-ii-renderer-red` and not on master, so it is that branch's to
correct, not this PR's. It is recorded here so that it does not propagate from
this file into `@tester`'s and `@developer`'s briefs, which is how the first
error spread.

**Gate A — the estimate gate, and it has already fired.** The measurement is the
per-file estimate in "Files and LOC" below, taken at design time. It comes to
~786 against `CLAUDE.md` §2's 700, so the split is **taken now** rather than
armed, and there is no threshold left to trip during 5b. **This is weaker
evidence than a measurement and should be read as such**: it is an estimate of
lines nobody has written, made by the person with the strongest interest in the
answer being "it fits". Its only virtue is that it is available at the one moment
the PR's shape is still free. Two things make it less weak than it sounds, and
neither makes it strong: it is built from per-file rows that a reader can
disagree with individually, and the two conversion factors applied to those rows
(C++ at or under, bindings at twice) are measured on four shipped increments
rather than guessed.

**Gate B — the shape gate, and it is `@tester`'s to run before the red
commit.** Not a line count, which is why it is the one gate here that can be
obeyed cheaply at the moment it fires:

> **One test file per production header. No test file may cover two headers.**
> Before committing the red suite, `@tester` checks the file list it is about to
> commit against the mapping in "What is worth testing" below, and raises a
> deviation **before** the commit.

The mapping is four production headers to four test files, plus the one declared
relaxation — `node.hpp`'s second, fixture file — which makes **five** registered
suites in all; see "What is worth testing". It is deliberately
one-to-one even where a shared file would be convenient — `broad_phase.hpp` and
`noded_pslg_builder.hpp` would sit naturally in one suite, and they are split
anyway. The cost of the discipline is a few extra `add_terrain_backend_test`
lines; the cost of not having it is the second red round 6b-ii declined to pay.

**Gate C — and it is a re-estimate instruction, not a trip-wire.** 5c's estimate
is ~341 and its largest row is `bindings/core.cpp` at ~215, the file whose
estimate has doubled twice. **Nothing measures that during 5c**, because by the
time there are lines to count the shape is fixed — so what is armed is an
obligation on whoever writes 5c's design: re-estimate `bindings/core.cpp` and
`_core.pyi` against what 5b's `NodeStatus`, `NodedPslg` and `NodeOutcome`
actually turned out to be, before `@tester` is briefed, and if the total passes
**550**, move `cli.py` and `--snap-spacing` to a 5d. Declared now so that 5c's
suite layout is chosen knowing it, and spelled as a duty rather than as a
threshold precisely because the threshold form is the one with the one-for-three
record.

## The types

### `core/noded_pslg.hpp` — the type, kernel-free

Namespace `terrain`. Includes `core/pslg.hpp` (for `Chain`, `ChainRole`,
`IndexedRing`) and `core/snap_grid.hpp`. It includes nothing from
`predicates/` and nothing from `noding/`: like `Pslg`, it is a type, and the
validation that establishes its promise is a separate, kernel-templated header,
exactly as `core/pslg.hpp` and `core/pslg_builder.hpp` are separate.

The surface is `05-noder.md`'s, fixed there and not re-litigated:

```cpp
// include/terrain/core/noded_pslg.hpp, namespace terrain
namespace noding { class NodedPslgBuilder; }

class NodedPslg {
public:
    [[nodiscard]] std::span<const Point2>        vertices() const noexcept;
    [[nodiscard]] std::span<const Chain>         chains() const noexcept;
    [[nodiscard]] std::span<const std::uint32_t> chain_indices() const noexcept;
    [[nodiscard]] std::span<const std::uint32_t> indices_of(std::size_t c) const noexcept;
    [[nodiscard]] IndexedRing                    ring(std::size_t c) const;
    [[nodiscard]] std::size_t                    edge_count(std::size_t c) const noexcept;
    [[nodiscard]] Segment2                       edge(std::size_t c, std::size_t k) const noexcept;

    [[nodiscard]] SnapGrid                       grid() const noexcept;
    [[nodiscard]] std::span<const EdgeProperties> edge_properties() const noexcept;
    [[nodiscard]] std::span<const std::uint32_t> node_of_input_vertex() const noexcept;
    [[nodiscard]] std::size_t                    edge_base(std::size_t c) const noexcept;

private:
    NodedPslg() = default;
    friend class noding::NodedPslgBuilder;
    // ... four vectors, one SnapGrid
};
```

Two additions to `05-noder.md`'s sketch, both forced:

- **`edge_base(c)`.** Guarantee 15 makes `edge_properties()` index-aligned with
  the flat enumeration `(c, k)`, and a caller holding `(c, k)` needs the offset
  to reach the byte. Without it every caller recomputes a prefix sum over
  `edge_count`, which is the kind of derived arithmetic that is right in four
  places and wrong in the fifth. `edge_base(0) == 0` and
  `edge_base(chains().size()) == edge_properties().size()` — one extra entry, the
  same convention `chain_indices` already uses for chain slices.
- **`friend class noding::NodedPslgBuilder`, not a friended function
  template.** `Pslg`'s precedent is `friend class PslgBuilder`
  (`include/terrain/core/pslg.hpp:164`), and friending `node<K>` instead would
  require `core/noded_pslg.hpp` to declare a `terrain::noding` function template
  parameterised on `pred::GeometryKernel` — dragging a concept from
  `predicates/` into a type header that has no business knowing kernels exist.

**`NodedPslg` is not a subclass of `Pslg` and there is no conversion**, per
`05-noder.md`. Restated because it is the whole architectural point: the
signature is what makes un-noded input unrepresentable at the CDT entry.

**`NodedPslg` satisfies `viz/protocols.py`'s `PslgLike` once bound** —
`vertices`, `chains`, `chain_indices`, `indices_of` — which is the renderer
requirement above and the reason `chains()` returns `Chain`, the same struct
`Pslg` uses, rather than a noded-specific one.

### Guarantees

On top of `Pslg`'s 1–8 and 10 (9 is replaced by `node_of_input_vertex`), and
renumbered from `05-noder.md`'s:

11. Every coordinate is `grid().world(g)` for some `GridPoint` with
    `|ix|, |iy| <= kMaxGridIndex`. Checkable bitwise: `grid().snapped(v) == v`.
12. No two vertices are equal.
13. No edge is zero-length; no chain has a repeated consecutive index.
14. **(a) No two edges cross in their interiors, and (b) no node's cell meets an
    edge it is not an endpoint of.** Amended from `05-noder.md`; see below.
15. `edge_properties()` is index-aligned with the flat edge enumeration
    `(c, k), k < edge_count(c)`, at offset `edge_base(c) + k`, and each entry is
    the **union** over the property sets of every input chain that contributed
    geometry to that edge. Union rather than a one-bit OR, per "Edge
    properties": an output edge can be a road *and* a river, and the merge is a
    reduce over an unordered set of contributors.
16. **Chain order and roles are preserved index-for-index with the input
    `Pslg`**: `chains().size()` equals the input's, and `chains()[c].role` and
    `chains()[c]`'s property set equal the input's. Only `begin` and `count`
    change.
    New here, and it is the mitigation for half of `05-noder.md` risk 5: per-
    *vertex* index identity is lost and `node_of_input_vertex` maps around it,
    but per-*chain* identity survives, so a Python-side parallel attribute array
    keyed by chain — which is what the chain property set itself is — stays valid without a
    map. A noder that dropped a collapsed ring silently would break this; it
    rejects instead (`RingCollapsed`).

#### Guarantee 14, amended

`05-noder.md` states 14(a) as "`classify<K>` returns `Disjoint` or `Touching`
for every candidate pair". **That is false on the road-along-a-river case, which
is the only input guarantee 15 does any work on.** After the split pass a road
and a river that run together occupy the same nodes, so both chains carry the
edge `(u, v)` with identical endpoints, and `classify` on two identical segments
returns `Overlapping`.

The amended clause:

> **14(a).** For every pair of distinct output edges, `classify<K>` returns
> `Disjoint` or `Touching`, **or** returns `Overlapping` and the two edges have
> **equal** node-id pairs `(min(u,v), max(u,v))`.

A partial overlap is still a violation, and after the split pass there is no
third case: if two collinear edges overlap at all, each one's endpoints lie in
the other's hot pixels, so each is split at the other's nodes and the overlap is
exact. Stating it as "equal node-id pairs" rather than "collinear" keeps the
check on integers, which is the shape `05-noder.md` rules the merge must keep.

14(b) is unchanged and is the hot-pixel clause: `segment_meets_cell<K>(grid,
edge, g)` is false for every (node, edge) pair with `g` not an endpoint of
`edge`. It is **strictly stronger than exact incidence**, which is why it
catches a T-junction the split pass missed instead of blessing it.

**Both clauses are established by verification, not by construction**, and the
verifier is a separate type so that it can be run against a broken input — see
`NodedPslgBuilder`.

#### The unit diagonal through a lattice corner: a worked degeneracy

`segment_meets_cell` uses **closed** cells and answers `true` on a grazed corner
(`include/terrain/noding/intersect.hpp:203-213`; 5a mutant 13 pins it). 5b
cannot change that, so 5b has to live with its consequence, which is sharper
than it looks:

Take nodes at grid indices `(0,0)`, `(1,1)` joined by an edge, and suppose
`(1,0)` and `(0,1)` are also nodes — some input vertex snapped to each. The edge
passes exactly through the lattice corner shared by all four cells. Both `(1,0)`
and `(0,1)` are therefore *grazed*, `segment_meets_cell` is true for both, and
guarantee 14(b) **requires** the edge to be split at both. Their dot products
along the edge are exactly equal, so the arc order is an exact tie, broken by
`GridPoint`'s defaulted `operator<=>`.

The output chain is `(0,0) → (0,1) → (1,0) → (1,1)` (or the mirror, depending on
the tie-break). Three things follow, and the third is the one that matters:

- It is **not** a doubling-back chain, so `05-noder.md` risk 10 does not apply
  here. The projections along the segment are `0, ½, ½, 1` — monotone
  non-decreasing. The zigzag is **transverse**, not longitudinal, and an exact
  tie can only ever be transverse: two nodes with equal projections differ
  perpendicular to the segment or not at all.
- The deformation is exactly `h/√2` — both inserted nodes sit on the boundary of
  the closed `h/√2`-neighbourhood — so the one-way Hausdorff bound
  (`05-noder.md`'s fix 8) holds, at equality. This is the fixture that shows the
  bound is tight and not merely an over-estimate.
- **Guarantee 14 does not hold on it. It was checked and it is false**, and the
  correction is the next subsection. An earlier revision of this bullet read
  "Guarantee 14 still holds on it, which `@tester` should check by hand rather
  than assume", and gave the reason: consecutive edges share an endpoint
  (`Touching`), and the two non-consecutive edges `(0,0)-(0,1)` and
  `(1,0)-(1,1)` are parallel and disjoint. Every clause of that is true and the
  conclusion is still wrong, because all of it is about **14(a)**, which quantifies
  over edge *pairs*. It never asks 14(b) — a (node, edge) claim — of the resolved
  chain's own interior nodes against its own middle edge. The bullet is kept in
  the record rather than deleted, because it is the shape of the error, not a typo:
  a clause that was checked against the wrong quantifier.

#### The tie does not converge, and the probe that shows it

`@tester` ran the check this section asked for and it came back red; it was then
re-run independently. The probe, against 5a's shipped
`segment_meets_cell<DefaultKernel>` and `SnapGrid{1.0}`, linking
`libterrain_predicates.a` — build directory per `CLAUDE.md` §4, and any of the
tree's `build*` directories carries the archive:

```cpp
// /tmp/probe5b.cpp
#include <cstdio>
#include "terrain/core/snap_grid.hpp"
#include "terrain/noding/intersect.hpp"
#include "terrain/predicates/default_kernel.hpp"
using namespace terrain;
using K = pred::DefaultKernel;
static void probe(const char* n, GridPoint a, GridPoint b, const SnapGrid& gr,
                  GridPoint c, GridPoint d) {
    Segment2 s{gr.world(a), gr.world(b)};
    std::printf("%s (%d,%d)-(%d,%d): meets (%d,%d)=%d  meets (%d,%d)=%d\n", n,
                (int)a.ix, (int)a.iy, (int)b.ix, (int)b.iy, (int)c.ix, (int)c.iy,
                (int)noding::segment_meets_cell<K>(gr, s, c),
                (int)d.ix, (int)d.iy,
                (int)noding::segment_meets_cell<K>(gr, s, d));
}
int main() {
    SnapGrid gr{1.0};
    probe("diag", {0,0}, {1,1}, gr, {1,0}, {0,1});   // the input edge
    probe("mid ", {0,1}, {1,0}, gr, {0,0}, {1,1});   // the resolved chain's middle edge
    probe("e1  ", {0,0}, {0,1}, gr, {1,0}, {1,1});   // its first edge
    probe("e3  ", {1,0}, {1,1}, gr, {0,0}, {0,1});   // its third edge
    return 0;
}
```

```sh
c++ -std=c++20 -Iinclude /tmp/probe5b.cpp build/libterrain_predicates.a -o /tmp/probe5b && /tmp/probe5b
```

prints

```
diag (0,0)-(1,1): meets (1,0)=1  meets (0,1)=1
mid  (0,1)-(1,0): meets (0,0)=1  meets (1,1)=1
e1   (0,0)-(0,1): meets (1,0)=0  meets (1,1)=0
e3   (1,0)-(1,1): meets (0,0)=0  meets (0,1)=0
```

Line 1 is why the edge is split at both. Line 2 is the finding: the resolved
chain's **middle edge is the anti-diagonal through the same lattice corner**, so
it grazes the closed cells of `(0,0)` and `(1,1)` — neither of which is an
endpoint of it. **That is a 14(b) violation in the output of round 1**, and
lines 3 and 4 say it is the only one: the two outer edges are clean.

So round 2 must split the middle edge at `(0,0)` and `(1,1)`. Their projections
along `(0,1) → (1,0)` are again exactly equal, the tie again resolves
lexicographically, and the middle edge becomes `(0,1), (0,0), (1,1), (1,0)` —
whose own middle edge is the **original diagonal**, which line 1 says grazes
`(1,0)` and `(0,1)`. Round 3 splits that. The two configurations map into each
other, the chain gains two positions per round and never satisfies 14(b):

```
round 0: (0,0) (1,1)
round 1: (0,0) (0,1) (1,0) (1,1)
round 2: (0,0) (0,1) (0,0) (1,1) (1,0) (1,1)
round 3: (0,0) (0,1) (0,0) (0,1) (1,0) (1,1) (1,0) (1,1)
```

**`Ok` is the one answer this input cannot receive**, at any `max_rounds`. The
terminating answer is `NotConverged`, and it is reached by the cap rather than by
a fixpoint. Note also that from round 2 the chain **doubles back** — `05-noder.md`
risk 10 after all, one round later than the bullet above says — and that the
repeated node pairs are duplicate edges, which amended 14(a) permits, so 14(a)
never fires and 14(b) is carrying the whole refusal.

A satisfying configuration does exist — lines 3 and 4 show that the L-shaped
chain `(0,0) → (0,1) → (1,1)` grazes nothing — and the split pass cannot reach
it, because the pass is defined to split an edge at **every** node whose closed
cell it meets, and 14(b) as written demands exactly that. The two are
inconsistent on an exact corner graze. This is a design defect, not a
`@tester` finding about an implementation, and it is ruled on as risk 16 below.

**Named fixture, `test_noding_node.cpp`: the lattice-corner tie.** Its purpose
has changed with the finding and is stated in full under risk 16; it is still
the cheapest input that exercises the tie-break, the grazed corner and the
tightness of the deformation bound at once, and it is now also the only fixture
that reaches `NotConverged`.

For a **closed ring** the same zigzag makes the ring non-simple by round 2 — the
repeated node ids above are precisely what the `NonSimpleRing` check looks for —
so a ring hits `NonSimpleRing` where an open chain hits `NotConverged`. Both are
risk 16; the input needed is a ring edge passing exactly through a lattice corner
flanked by two nodes, at one cell's separation.

### `noding/broad_phase.hpp` — the index, and its one query

Namespace `terrain::noding`. Includes `core/bbox.hpp` and `core/segment.hpp`.
**It does not include `core/snap_grid.hpp`**, per `05-noder.md`, and that is
grep-checkable.

```cpp
// Visits every indexed segment whose closed bbox intersects q, and may visit
// others. F is invoked as f(std::uint32_t segment_index).
template <class F> void for_each_candidate(Box2 q, F&& f) const;
```

The contract is `05-noder.md`'s and is not weakened: *the visited set contains
every segment whose closed bounding box intersects `q`*. False positives are
free; a false negative is a bug in a function whose entire job is not to have
them. **No `Point2` overload exists**, for the reason `05-noder.md` gives — the
one live way to break conservativeness is querying a node by its centre instead
of by its cell's box.

Two things 5b fixes that `05-noder.md` left open:

**Bucketing is by bounding box, not by exact traversal, and that is what makes
the postcondition provable in one line.** Each segment is inserted into every
bucket its closed bbox overlaps. A query box `q` that intersects segment `s`'s
bbox `B` shares at least one bucket with `B`, because both are rasterised by the
same bbox-overlap rule over the same lattice — so `s` is visited. An exact
traversal (a DDA along the segment) would be *tighter* and would make that
argument a case analysis instead of a sentence. `05-noder.md` forbids bucketing
by endpoints; it does not require exactness, and exactness here buys fewer false
positives at the cost of the only correctness argument that fits on a line.

**Bucket sizing is by segment count and domain extent, never by
`SnapGrid::spacing()`.** `k = max(1, ceil(sqrt(n)))` buckets per axis over the
segment set's bounding box, so occupancy is O(1) for uniformly distributed
segments and the index is O(n) in memory. Degenerate extents are handled without
a branch in the query: a zero-width or zero-height extent gives one bucket on
that axis. Non-finite coordinates cannot arrive — `Pslg` guarantee 1.

The index is **immutable after construction** and const-queried, which is what
lets the split pass and the verification pass share one instance and lets a
future parallel-for over segments query it without synchronisation.

### `noding/noded_pslg_builder.hpp` — the verifier, and the only producer

Namespace `terrain::noding`. This header is **not** in `05-noder.md`'s shipping
list; it is the contingency split's `NodedPslgBuilder` promoted to a first-class
part of 5b, and the reason is a testability one that outranks the line count:

> **A verification pass buried inside `node<K>`'s loop cannot be run against a
> broken input.** The only way to reach it would be through the noder, so the
> only inputs it ever sees are the ones the noder produced, and a verifier that
> agrees with its producer is `.claude/REQUIRED-READING.md`'s "do not write X is
> verified by Y until Y has been run against a broken X" — the exact defect
> `05-noder.md` guarantees 14 and 15 were revised three times to remove.

As a separate type with a public entry point, the verifier can be handed a
hand-built candidate with a live T-junction in it and be *required to refuse*.
That is the only suite in this increment that can fail when the verifier is
broken, and it is why the builder is invariant-critical.

```cpp
// include/terrain/noding/noded_pslg_builder.hpp, namespace terrain::noding
enum class NodeStatus : int { Ok = 0, NotRun, InvalidSnapSpacing, CoordinateOutOfRange,
                              RingCollapsed, RingDegenerateAfterSnap, NonSimpleRing,
                              NotConverged, MalformedOutput };

[[nodiscard]] constexpr std::string_view describe(NodeStatus) noexcept;

struct NodeOutcome {
    std::optional<NodedPslg> pslg;   // engaged iff status == Ok
    NodeStatus  status{NodeStatus::NotRun};
    std::string message;

    [[nodiscard]] bool ok() const noexcept { return pslg.has_value(); }
};

class NodedPslgBuilder {
public:
    // Takes the candidate arrays by value and consumes them. The kernel is a
    // parameter of build(), not of the class, exactly as PslgBuilder::build is.
    NodedPslgBuilder(SnapGrid, std::vector<Point2>, std::vector<Chain>,
                     std::vector<std::uint32_t>, std::vector<EdgeProperties>,
                     std::vector<std::uint32_t>);

    template <pred::GeometryKernel K>
    [[nodiscard]] NodeOutcome build() &&;
};
```

`build()` checks guarantees 11, 12, 13, 14(a), 14(b) and 16 and produces a
`NodedPslg` only if all hold. It does **not** check 15 — see below.

**The minimum-node-count guard is `is_closed(role)`'s alone**, so an open chain
of a single node is accepted and an open chain is never asked for three. That is a
ruling rather than a consequence of how the condition was written; it is argued
and pinned under "An open chain that collapses to one node is `Ok`".

**Ordering, and it is not arbitrary.** 11 and 12 come first because they are
O(n) and because 14's node-id reasoning is meaningless without them; 13 next; 14
last, because it is the expensive one. 16 is a comparison against the input
chain table and so is the caller's to supply — the builder receives the chain
vector and checks internal consistency (contiguous slices, no gaps), and the
*correspondence* to the input is `node<K>`'s to establish and the unit suite's
to check.

**`MalformedOutput` is a self-check with a name**, exactly parallel to
`CdtStatus::MalformedInput` (`include/terrain/cdt/result.hpp:55-59`). Guarantees
11, 12, 13 and 16 are things `node<K>` establishes by construction; the builder
refusing one of them means the driver is broken, not that the input was bad.
14 is different: it is *not* established by construction, snap rounding can
genuinely create a crossing that was not in the input, and the driver's response
is to iterate. So the driver maps a 14 refusal to `NotConverged` when the cap is
reached and to another round otherwise, and maps 11/12/13/16 straight through as
`MalformedOutput`. **Two statuses for one refusal, because they mean opposite
things to the person reading them**: `NotConverged` says "try a finer spacing",
`MalformedOutput` says "this is our bug".

**And `build()` itself must name them, which an earlier revision left to the
driver alone.** `@tester` found the gap: the driver cannot map what it cannot
distinguish, and `build<K>()` returns a single `NodeStatus`, so the discriminator
has to be *in* that status or nowhere. Confirmed as `@tester` pinned it:

> **`NodedPslgBuilder::build<K>()` returns `NodeStatus::NotConverged` for a
> refusal of guarantee 14 (either clause) and `NodeStatus::MalformedOutput` for a
> refusal of 11, 12, 13 or 16.**

The reading that makes that coherent, and it belongs at the head of the header
rather than in a suite comment: **`NodeStatus` is the status of an *outcome*, not
of a check. The builder reports the status the driver would report if this were
the last round.** Under that reading `build()`'s mapping is the identity at the
cap — the driver forwards rather than translates — and the only place the name
reads oddly is the unit suite, where a hand-built T-junction candidate comes back
`NotConverged` having converged nothing. That oddity is confined to one file read
by people who have this paragraph; a third enumerator to remove it would cost an
enum row, a `describe` arm, a Python enum member at 5c and a mapping the driver
would immediately collapse. The nine-row enum stands.

A refusal that violates **both** 14(b) and 11/12/13/16 reports `MalformedOutput`:
the self-check outranks the diagnosis, because a builder that is being handed
malformed arrays has no basis for saying anything about convergence. 14(a) and
14(b) both reporting is not such a case — both are 14, both are `NotConverged`,
and see the partial-overlap note under the builder's suite.

**Why guarantee 15 is not the builder's.** Its oracle must be built from the
*input* — that is `05-noder.md`'s guarantee-15 ruling, arrived at after the
obvious check turned out to be the dedup restating itself — and the builder does
not have the input. Verifying 15 inside the builder could only mean re-deriving
it from the arrays the driver just handed over, which is the self-confirming
shape by another route. **15 is checked in the property suite, against the
input, and nowhere else.** That is a real, stated gap in what the *type* system
guarantees, and the honest form of it: a `NodedPslg` promises 11–14 and 16 by
construction, and 15 only as far as the suite reaches.

### `noding/node.hpp` — the driver

Namespace `terrain::noding`. Includes all of the above plus
`noding/intersect.hpp`, `noding/node_set.hpp`, `core/pslg.hpp` and
`core/ring.hpp`.

```cpp
struct NodeOptions {
    double        spacing;            // no default, ever; see 05-noder.md risk 2
    std::uint32_t max_rounds{4};
};

template <pred::GeometryKernel K>
[[nodiscard]] NodeOutcome node(const Pslg& pslg, const NodeOptions& options);
```

`spacing` has **no default** and is the first member so that
`NodeOptions{0.1}` is the only cheap spelling. `05-noder.md` risk 2 rules that
there is no default spacing anywhere in C++, for the same reason `K` has no
default in `ring.hpp` and `pslg_builder.hpp`. A default in `cli.py` at 5c is a
policy choice by the composition root and is a different thing.

`max_rounds{4}` **does** have a default, because it is a cap on a loop rather
than a parameter of the answer: at any value the outcome is either the same mesh
or `NotConverged`, never a different mesh. Four rather than one because
`05-noder.md` risk 1 is real, and four rather than sixteen because the audit
measured the cascade to be rarer than the theory allows (`05-noder.md` risk 9:
the near-parallel generator produced no surviving crossings at 0.1 m). **The
design does not claim convergence in four rounds**; `NotConverged` is a
legitimate outcome with an actionable lever behind it.

`node<K>` is a **pure function** of `(pslg, options)`. No statics, no caches, no
global state, no allocation outside its own locals. Two calls on the same input
from two threads produce bit-identical output, because every ordering in it is
either `GridPoint`'s total order or the input's chain order, and never the order
a thread found something in. That is what makes it safe to release the GIL
around at 5c and what `05-noder.md` mutant 12 — an input-dependent anchor — would
destroy.

#### The split pass

`05-noder.md`'s six steps, unchanged in substance, with the things 5b has to
decide made explicit:

0. **Admission.** `is_valid_spacing(options.spacing)` → `InvalidSnapSpacing`.
   One pass over `pslg.vertices()` with `SnapGrid::can_snap` →
   `CoordinateOutOfRange`, naming the vertex. This is the one place
   `kMaxGridIndex` is compared against, per `snap_grid.hpp`'s comment, and after
   it `snap` is branch-free.
1. **Broad phase** over the current segment set → candidate segment pairs. Each
   segment's bbox is the query.
2. **`classify<K>` per pair.** `Crossing` → `crossing_point<K>` contributes a new
   node. `Touching` and `Overlapping` contribute nothing new. This is the only
   step that manufactures a coordinate.
3. **`NodeSet`** over every snapped input vertex and every constructed crossing.
   Node ids are its sorted order and are therefore a function of the *set*.
4. **Broad phase again, segments against nodes**, queried by the node's
   `cell_min`/`cell_max` box — never by `world(g)`. For each segment the nodes
   `g` with `segment_meets_cell<K>(grid, seg, g)`, sorted by the dot product of
   `(world(g) - s.a)` with `(s.b - s.a)`, ties broken by `GridPoint`'s `<=>`, are
   its split points. **The endpoints are anchored, not sorted** — see below.

##### The split sequence is anchored at its own endpoints

An earlier revision of step 4 said only "sorted by the dot product", and left the
implementation to discover that a plain sort of the met set is not a split
sequence. It is ruled here because the choice is the design's:

> **A segment's split sequence is `node(s.a)`, then the nodes it meets in arc
> order with `node(s.a)` and `node(s.b)` removed, then `node(s.b)`.** The
> endpoints are placed, never sorted into position. A segment whose two endpoints
> snap to the same node has the one-element sequence `[node(s.a)]` and
> contributes no edge.

The reason is continuity, and it is a consequence of the cell's reach rather than
of anything the sort does wrong. A cell extends `h/√2` from its centre, so a node
`x` grazed near `s.a` can have `dot(world(x) - s.a, s.b - s.a) < 0` — it projects
*behind* the segment's own start. A plain sort then puts `x` first, the segment's
piece list begins `x → node(s.a)`, and the reassembled chain no longer begins at
the node the previous segment ended at. The chain is discontinuous and nothing
downstream says so: guarantees 11–14 are all satisfied by a discontinuous run,
and `NonSimpleRing` only sees repeated ids. The same argument holds at `s.b`.
Anchoring costs one `push_back` at each end and removes the failure mode
entirely. `node.hpp` carries the rule as a comment where the sequence is built.

**Does it agree with the plain sort?** In every non-degenerate case, yes, and the
two spellings are then textually different and observationally identical:
`node(s.a)` and `node(s.b)` are always *members* of the met set — each endpoint
lies on the segment and inside its own node's cell, so `segment_meets_cell`
admits it unconditionally — and when no grazed node projects outside `[0, |s|]`
they are already the arc-order extremes. The spellings diverge on exactly the
inputs where one does: the anchored sequence keeps `x` interior, the plain sort
promotes it to an end.

**And the divergence is observable, which is why the oracle is not left alone.**
With met arc-sorted as `[x, a, y, b]` (`a = node(s.a)`, `b = node(s.b)`), the
anchored sequence is `[a, x, y, b]` and yields the undirected pairs `{a,x}`,
`{x,y}`, `{y,b}`; the plain sort yields `{x,a}`, `{a,y}`, `{y,b}`. `{a,x}` and
`{x,a}` are the same pair, so the symmetric difference is **two** pairs, not one:
the plain sort **gains** `{a,y}`, which the driver never produces, and **loses**
`{x,y}`, which it does. An earlier revision of this paragraph said "exactly one
pair" and was counting only the half it went on to argue from.

**The ruling survives because the two directions are not symmetric, and that
asymmetry is the whole argument.** `check_guarantee_15` only ever *unions* into
`proven` and then asserts `out.edge_properties()[e.flat].contains(proven)`. So:

- a **wrongly-claimed** pair can push `proven` past what the edge's properties
  really are — if some other chain produced `{a,y}`, this chain's properties are
  unioned in for an edge it never touched, `proven` stops being a subset, and the
  `REQUIRE` goes **red on correct output**. That is increment 5's exact-incidence
  failure with a different cause;
- a **wrongly-omitted** pair merely makes `proven` smaller, which an assertion
  that is already one-directional absorbs by construction. It is permanently
  inert, not inert-for-now.

So the correct claim is **exactly one pair in the direction that can fail**, and
it is the gained one. The lost pair is why the oracle is weaker than it looks,
not why it is wrong.

**And it is latent in the strongest sense: the generator never produces the
case.** `@tester` instrumented `split_sequence` to compute both spellings side by
side across the suite's generated sets and found **0 divergences over ~1000 input
segments**. So the grazed-node-behind-`s.a` case is not merely inert under the old
oracle — no seed reaches it, and no amount of running the old suite would have
found this. It was found by reading the rule against the code, which is the
`is the claim about the same object the code evaluates?` inspection rather than a
run, and is the case `.claude/REQUIRED-READING.md` reserves that question for.

> **Ruling: `split_sequence` in `tests/cpp/property/prop_noding_no_crossings.cpp`
> should be respelled to anchor its endpoints, matching this rule.** That is
> borrowing the producer's *specification*, not its records — the rule above is
> normative here, so an oracle written to it still refutes a driver that deviates,
> whereas the plain sort encodes a spec this design does not state. It is
> `@tester`'s file and lands as `@tester`'s own commit with the reason in the
> message, per `docs/increments/README.md`. `@architect` does not touch it.
5. **Reassemble chains and merge.** Each input chain becomes one output chain
   with the same role and property set (guarantee 16), its index run being the
   concatenation of its segments' node sequences with the shared endpoints
   written once. Zero-length edges — two consecutive positions with the same
   node id — collapse. **The chain structure is preserved; the dedup is over the
   edge *set*, and what it produces is the property set.** For each flat
   position `(c, k)` the entry is the union over every input chain carrying an
   edge with the same node-id pair.
6. **Verify** by handing the arrays to `NodedPslgBuilder::build<K>()`. `Ok` →
   done. A 14 refusal below the cap → iterate from 1 on the *split* segments. A
   14 refusal at the cap → `NotConverged`. An 11/12/13/16 refusal →
   `MalformedOutput`, immediately, without another round.

**"No new node" is not the fixpoint and must not be the loop condition.** A
round can split an edge at an existing node — a new edge, no new node — so a
loop spelled that way reports converged in a state that violates 14(b).
`05-noder.md` fixes this and it is restated here because it is the single
likeliest way for this increment to ship a surviving T-junction. It is mutant 5
below.

#### The step that is not in `05-noder.md`'s list: ring re-validation

Between 5 and 6, for every closed chain:

- **Distinct node count.** Fewer than 3 → `RingCollapsed`, naming the chain and
  the spacing.
- **Winding under `K`**, re-evaluated on the *snapped* ring with
  `orientation<K>` (`include/terrain/core/ring.hpp:233`). An `Outer` ring that is
  no longer counterclockwise, or a `Hole` no longer clockwise, or either now
  `Collinear` → `RingDegenerateAfterSnap`. **Never silently reversed**, following
  increment 3's stage 5 rule. **It is a self-check, not a diagnosis** — see
  "`RingDegenerateAfterSnap` is a self-check" below, which is the ruling and also
  relocates `05-noder.md` risk 3's mitigation.
- **Simplicity**, in the one form the noder can decide for free: a closed chain
  that visits the same node id twice. → `NonSimpleRing`. This is the
  figure-eight case, and it is detected by a `std::vector<std::uint32_t>` sort
  rather than by geometry, because after noding a ring self-intersects **iff** it
  repeats a node — 14(a) has already excluded interior crossings that are not at
  nodes. That equivalence is what makes the check cheap and it is worth a
  comment in the header, because the obvious implementation is a quadratic
  segment-pair scan.

#### `RingDegenerateAfterSnap` is a self-check, and risk 3's mitigation moves

`@tester` measured four candidate slivers designed to flip a hole's winding under
snapping. **All four were refused as `NonSimpleRing` before the winding check
could fire**, and reordering the throwaway's checks into this section's stated
order (count → winding → simplicity) did not change it, because the refusal is
not raised by the ordering — it is raised one stage earlier, by the split pass.

The mechanism, and it is exact rather than statistical. Take a ring whose
snapped winding has flipped or collapsed. Its snapped signed area has crossed
zero, which puts some vertex `C` on the far side of the chord its neighbours span
by less than the rounding scale — so `C` lies within half a cell of an edge it is
not an endpoint of, which is precisely `segment_meets_cell`. Guarantee 14(b) then
*requires* that edge to be split at `C`, and the ring `… A, C, B …` with `A–B`
split at `C` becomes `… A, C, B, C …`. **It repeats a node id, which is the
`NonSimpleRing` predicate verbatim.** The collinear case is the same argument with
no slack at all: any ring whose nodes are exactly collinear has a spanning edge
containing the others, is split at them, and repeats them.

**Ruling: `RingDegenerateAfterSnap` is demoted to a self-check with a name,
exactly like `MalformedOutput` and `CdtStatus::NotNoded`.** It stays in the
enum and the check stays in the code — deleting a cheap check whose job is to be
unreachable is how a later change makes it reachable in silence — but its
`describe` text says "self-check", its status-table lever column says so, and
**no fixture is owed for it**, because none can exist. That last clause is the
evidence for the demotion rather than an excuse for the gap: a status with no
reachable input is a status with no test, and the honest form is to say which
kind of status it is.

#### The ruling depends on an increment 2 implementation choice, and nothing says so

The paragraph above says reordering the checks into count → winding → simplicity
did not change the outcome. That holds only because the winding check does not
fire on a ring the split pass has just made non-simple — and **whether it fires
depends on which of two instruments `ring.hpp` spells `orientation<K>` as**. The
repo contains both, and on exactly the shape a split produces they disagree.

Probed against the shipped headers, compiled and run rather than reasoned about.
Take `A = (0,0)`, `B = (10,0)`, `C = (5,1)` — the ring `A,B,C` with `A–B` split at
a grazing `C`, which is the `… A, C, B, C …` shape the mechanism above
manufactures:

```
orientation<FastKernel>(A,B,C,B) = Clockwise     signed_area(A,B,C,B) = 0
orientation<FastKernel>(A,C,B,C) = CounterClockwise  signed_area(A,C,B,C) = 0
```

The zero is not a rounding artefact and does not depend on the coordinates: with
the shoelace translated to `vertex(0)`, the ring `A,X,B,X` contributes exactly
`cross(X-A, B-A)` and `cross(B-A, X-A)`, which are exact negatives, so the sum is
identically zero for **every** split sliver of this shape. A shoelace sign would
return `Collinear`, `Collinear != want`, and the ring would raise
`RingDegenerateAfterSnap` — the status this increment has just demoted as
unreachable — on the very input the demotion is argued from.

`orientation<K>` survives because increment 2 spelled it as an **extreme-vertex
walk advancing `prev` and `fwd` independently** (`include/terrain/core/ring.hpp:233`,
with the reason in the comment above it: a lockstep walk returns `Collinear` for
the proper triangle `{A,A,B,C}`). A repeated vertex costs it one collinear triple
and it walks past. Measured, because "walks past" is a claim: over 500 000
randomly generated rings of the shape a split produces — vertex `D` placed within
0.05 of edge `A–B` and `A–B` split at it — and separately over 3 200 000 rings
with a repeated vertex inserted at an arbitrary position, `orientation<FastKernel>`
returned `Collinear` **zero times**. `signed_area` is zero on all of the
sliver-shaped ones by the identity above.

**What would break it, stated plainly: anyone "simplifying" `orientation<K>` to a
shoelace sign silently reopens a status this increment closed.** The two look
interchangeable — both answer "which way does this ring wind" — and the shoelace
is shorter, cheaper and kernel-free, which is the whole reason the simplification
is tempting. `ring.hpp`'s own comment already forbids it for a different reason
(cancellation over a sliver at UTM33 magnitudes); this is a second, independent
reason, and it lives here because this is where the ruling that depends on it is.
`docs/increments/02-core-geometry.md` may deserve the same note beside
`orientation`'s own design — **that is another file and this section does not edit
it**; whoever next opens it can carry the probe across.

**One residual, named rather than absorbed.** `orientation<K>` returning a
*definite* orientation is measured above and is what the demotion needs. Whether
it returns the orientation matching the chain's `role` is a separate question and
is **not** invariant: in the 500 000-ring probe the split ring's reading differed
from the pre-split ring's in 3 410 cases. On those shapes the winding check would
fire first and report `RingDegenerateAfterSnap` where the honest diagnosis is
`NonSimpleRing`. That does not make the status reachable — the input must also be
a valid `Pslg` whose *snapped* ring takes such a shape, and the bounded search
below found no surviving flip at all — but it is the clause the four measured
slivers do not cover, and the "no fixture can exist" claim above is therefore
*not refuted by the search*, not *proved*.

**It is 5d's to settle, not this increment's, and that is a ruling rather than a
deferral.** At 5b the 3 410 in 500 000 changes *which name a refusal carries* —
`RingDegenerateAfterSnap` where `NonSimpleRing` is the honest diagnosis — and not
whether the ring is refused or what the caller should do about it. Nobody is owed
that distinction here. At 5d they are, because 5d narrows
`segment_meets_cell<K>`, which is the predicate the whole demotion rests on: the
argument is "the split pass refuses it one stage earlier", and what decides
whether the split pass fires is exactly that predicate. `@tester` and
`@orchestrator` reached this independently and it is recorded in
`docs/increments/05d-corner-graze.md` under "`RingDegenerateAfterSnap`'s
demotion rests on this predicate", **with two obligations rather than one**: the
one-line orientation report *and* a re-run of the whole bounded ring search
against the new predicate. The second is the one at risk of being skipped —
`SURVIVORS 0` is a measurement of the *old* predicate, so re-running only the new
probe tests the old question with the new code.

**And `05-noder.md` risk 3's mitigation is not where that document thinks it
is.** Risk 3 reads: if the winding re-check is dropped, "a reversed sliver hole
reaches `detria` and meshes as an island with no diagnostic anywhere". That is
false on this design. Drop the winding re-check entirely and the same ring is
still refused — by `NonSimpleRing`, from the node repeat the split pass
manufactured, one stage earlier. The mitigation for risk 3 is **guarantee 14(b)
plus the node-repeat check**, and the winding re-check is defence in depth behind
it. `05-noder.md`'s entry is corrected accordingly in this PR; the risk itself is
real and unchanged, only its mitigation moves.

**The collinear and sliver cases are settled by argument; the residual was
settled by search.** The collinear case needs no code — a ring whose nodes are
exactly collinear has a spanning edge containing the others. The sliver case has
the area bound behind it: a simple ring's signed area changes by at most its
perimeter times `h/√2` under snapping, so a flip needs area below that, which is
the sliver regime. What neither settles is the residual — **that no ring exists
which is simple before snapping, simple *and* node-repeat-free after it, and
reversed**, since such a ring would reach the winding check and make
`RingDegenerateAfterSnap` a diagnosis after all. That one is not an argument's to
settle, so it was searched.

#### The bounded ring search, and the intermediate result that looks like a refutation

Read each block's stages in order, because **the flip count alone looks like a
counterexample set and is not one**. A search that stops at "the winding
flipped" reports thousands of refutations of the ruling and exits non-zero.
Applying the mechanism's own filter — does any node's cell meet an edge it is not
an endpoint of — takes every one of them to zero. **`FLIPPED` equals `grazed`
and `SURVIVORS` is zero in every block**, which is the claim, as a number rather
than as a structural argument:

```
4-gons, 0.25 step, 9x9, loose simplicity, simple-after REQUIRED:
  tested 2459896 | repeat-id 1376569 | not-simple-after 301222 | same-winding 781252 | FLIPPED 853 | grazed 853 | SURVIVORS 0
4-gons, 0.25 step, 9x9, loose simplicity, simple-after not required:
  tested 2459896 | repeat-id 1376569 | not-simple-after 0 | same-winding 1080766 | FLIPPED 2561 | grazed 2561 | SURVIVORS 0
4-gons, 0.25 step, 9x9, strict simplicity, simple-after REQUIRED:
  tested 2459896 | repeat-id 1376569 | not-simple-after 301222 | same-winding 781252 | FLIPPED 853 | grazed 853 | SURVIVORS 0
5-gons, 0.50 step, 6x6, loose simplicity, simple-after REQUIRED:
  tested 1217716 | repeat-id 538504 | not-simple-after 135336 | same-winding 541264 | FLIPPED 2612 | grazed 2612 | SURVIVORS 0
5-gons, 0.50 step, 6x6, loose simplicity, simple-after not required:
  tested 1217716 | repeat-id 538504 | not-simple-after 0 | same-winding 673192 | FLIPPED 6020 | grazed 6020 | SURVIVORS 0
5-gons, 0.50 step, 6x6, strict simplicity, simple-after REQUIRED:
  tested 1016516 | repeat-id 435812 | not-simple-after 131636 | same-winding 449068 | FLIPPED 0 | grazed 0 | SURVIVORS 0
5-gons, 0.25 step, 7x7, loose simplicity, simple-after REQUIRED:
  tested 6263524 | repeat-id 5457883 | not-simple-after 191199 | same-winding 614111 | FLIPPED 331 | grazed 331 | SURVIVORS 0
```

The four stages plus `FLIPPED` sum to `tested` in every block, which is what says
no ring is silently dropped: a rejection stage that quietly ate candidates would
show up as a shortfall, and that arithmetic is the only self-check the program
has. The program, inline — a design that names a throwaway it does not contain
is not reproducible by the next reader, whoever still has a copy today:

```cpp
// /tmp/ringsearch.cpp
// Bounded exhaustive search: does a ring exist that is simple and non-degenerate
// BEFORE snapping, has no repeated snapped node id and is still simple AFTER,
// and whose winding has flipped -- WITHOUT any node's cell meeting an edge it is
// not an endpoint of? Such a ring would be reachable RingDegenerateAfterSnap.
#include <algorithm>
#include <cstdio>
#include <vector>
#include "terrain/core/snap_grid.hpp"
#include "terrain/noding/intersect.hpp"
#include "terrain/predicates/default_kernel.hpp"

using namespace terrain;
using K = pred::DefaultKernel;
using noding::SegmentRelation;

static double shoelace(const std::vector<Point2>& r) {
    double t = 0.0;
    const Point2 o = r[0];
    for (std::size_t i = 1; i + 1 < r.size(); ++i) t += cross(r[i] - o, r[i + 1] - o);
    return 0.5 * t;
}

// Simple: adjacent edges meet only at their shared endpoint (Touching, never
// Overlapping), non-adjacent edges are Disjoint. Exact, via 5a's classify<K>.
static bool simple(const std::vector<Point2>& r, bool strict) {
    const std::size_t n = r.size();
    for (std::size_t i = 0; i < n; ++i) {
        const Segment2 a{r[i], r[(i + 1) % n]};
        if (a.a.x == a.b.x && a.a.y == a.b.y) return false;
        for (std::size_t j = i + 1; j < n; ++j) {
            const Segment2 b{r[j], r[(j + 1) % n]};
            const bool adj = (j == i + 1) || (i == 0 && j == n - 1);
            const SegmentRelation rel = noding::classify<K>(a, b);
            if (adj) {
                if (rel != SegmentRelation::Touching) return false;
            } else if (strict) {
                if (rel != SegmentRelation::Disjoint) return false;
            } else {
                // Loose: proper crossings only. A vertex merely touching a
                // non-adjacent edge is admitted.
                if (rel == SegmentRelation::Crossing) return false;
            }
        }
    }
    return true;
}

// The mechanism under test: 14(b) over the ring's own nodes and its own edges.
static bool grazes(const SnapGrid& g, const std::vector<GridPoint>& ids,
                   const std::vector<Point2>& r) {
    const std::size_t n = r.size();
    for (std::size_t e = 0; e < n; ++e) {
        const Segment2 s{r[e], r[(e + 1) % n]};
        for (std::size_t v = 0; v < n; ++v) {
            if (v == e || v == (e + 1) % n) continue;
            if (noding::segment_meets_cell<K>(g, s, ids[v])) return true;
        }
    }
    return false;
}

struct Counts { long tested = 0, repeat = 0, nonsimple_after = 0, same = 0,
                flipped = 0, grazed = 0, survivors = 0; };

static void run(const char* label, std::size_t n, double step, int side, bool strict,
                bool require_simple_after) {
    const SnapGrid g{1.0};
    std::vector<Point2> cand;
    for (int i = 0; i < side; ++i)
        for (int j = 0; j < side; ++j) cand.push_back(Point2{i * step, j * step});

    Counts c;
    std::vector<std::size_t> pick(n);
    const std::size_t m = cand.size();
    // Every n-subset; every distinct cyclic order of it (first index pinned,
    // reflections halved). Recursion-free combination walk.
    for (std::size_t i = 0; i < n; ++i) pick[i] = i;
    for (;;) {
        std::vector<std::size_t> rest(pick.begin() + 1, pick.end());
        std::sort(rest.begin(), rest.end());
        do {
            if (rest.back() < rest.front()) continue;  // drop the reflection
            std::vector<Point2> pre{cand[pick[0]]};
            for (std::size_t k : rest) pre.push_back(cand[k]);
            if (!simple(pre, strict)) continue;
            const double a0 = shoelace(pre);
            if (a0 == 0.0) continue;

            std::vector<GridPoint> ids;
            std::vector<Point2> post;
            for (const Point2& p : pre) {
                ids.push_back(g.snap(p));
                post.push_back(g.world(ids.back()));
            }
            bool repeat = false;
            for (std::size_t x = 0; x < n && !repeat; ++x)
                for (std::size_t y = x + 1; y < n && !repeat; ++y)
                    repeat = (ids[x] == ids[y]);
            ++c.tested;
            if (repeat) { ++c.repeat; continue; }   // refused earlier: repeated node id
            if (require_simple_after && !simple(post, strict)) {
                ++c.nonsimple_after; continue;     // refused earlier: crossing route
            }
            const double a1 = shoelace(post);
            if (a1 == 0.0 || (a1 > 0.0) == (a0 > 0.0)) { ++c.same; continue; }
            ++c.flipped;
            if (grazes(g, ids, post)) ++c.grazed; else ++c.survivors;
        } while (std::next_permutation(rest.begin(), rest.end()));

        std::size_t i = n;
        while (i-- > 0) if (pick[i] != i + m - n) break;
        if (i == static_cast<std::size_t>(-1)) break;
        ++pick[i];
        for (std::size_t j = i + 1; j < n; ++j) pick[j] = pick[j - 1] + 1;
    }
    std::printf("%s, %.2f step, %dx%d, %s simplicity, simple-after %s:\n"
                "  tested %ld | repeat-id %ld | not-simple-after %ld | same-winding %ld"
                " | FLIPPED %ld | grazed %ld | SURVIVORS %ld\n",
                label, step, side, side, strict ? "strict" : "loose",
                require_simple_after ? "REQUIRED" : "not required",
                c.tested, c.repeat, c.nonsimple_after, c.same, c.flipped, c.grazed,
                c.survivors);
}

int main() {
    run("4-gons", 4, 0.25, 9, false, true);
    run("4-gons", 4, 0.25, 9, false, false);
    run("4-gons", 4, 0.25, 9, true, true);
    run("5-gons", 5, 0.50, 6, false, true);
    run("5-gons", 5, 0.50, 6, false, false);
    run("5-gons", 5, 0.50, 6, true, true);
    run("5-gons", 5, 0.25, 7, false, true);
    return 0;
}
```

```sh
c++ -std=c++20 -O2 -Iinclude /tmp/ringsearch.cpp build/libterrain_predicates.a -o /tmp/ringsearch && /tmp/ringsearch
```

Simplicity is decided by 5a's own `classify<K>` rather than by a hand-rolled
segment test, so the search borrows the *predicate* and none of the producer's
records, per `computational-geometry/SKILL.md`.

**Two questions, two sets, and the design picks neither — it labels both.** The
`simple-after` switch is not a tuning knob, it selects which question is being
asked, and a reader who reproduces one block and compares it to the other's
number will think they have found something:

- **`simple-after REQUIRED` is the set of candidate counterexamples**, and it is
  the right set for "is `RingDegenerateAfterSnap` reachable". A ring that is
  non-simple after snapping is already refused by `NonSimpleRing` through the
  *crossing* route, so it never needed the graze argument at all and cannot be a
  counterexample to it.
- **`simple-after not required` is the stronger test of the graze mechanism
  itself**, precisely because it admits rings the tight set excludes — 6,020
  five-gon flips against 2,612 — and every one of them still grazes. That is a
  harder test of the mechanism, not a larger count of counterexample candidates,
  and the two claims must not be swapped.

**Two other instrument findings, both of which cut against the obvious reading.**

*The strict simplicity filter measures nothing on 5-gons.* Requiring non-adjacent
edges to be `Disjoint` excludes rings where a vertex merely touches a
non-adjacent edge, and it is so tight that the 5-gon block admits **zero** flips
— output indistinguishable from a search that measured everything and found
nothing, which is exactly the hazard `.claude/REQUIRED-READING.md` names. The
loose filter, rejecting proper crossings only, produces the 2,612. On 4-gons the
two agree exactly, which is the check that the loose one has not been loosened
into admitting garbage.

*A finer sub-cell step makes flips rarer, not commoner.* The 5-gon block at a
0.25 step over 7×7 tests 6,263,524 rings — five times the 0.50 block — and finds
**331** flips against 2,612, because finer steps put more vertices in one cell
and `repeat-id` rejects 5,457,883 of them, 87 % of the total, before the winding
is ever compared. **The flip count is not a measure of how hard the search
tried.** A reader tuning the lattice to maximise flips is tuning toward the wide
set, not toward a better search.

**Two independently written programs, and what is and is not reconciled.** This
program is a reconstruction from a stated method, written without sight of the
original; the original was re-run with the same stage instrumentation. The
load-bearing result is identical in both, and so is the 87 % `repeat-id`
rejection at the finer lattice, which is a specific enough agreement to say the
two programs are measuring the same thing. **The totals are not reconciled.**
The `simple-after` filter is a real part of the gap and it is not the whole of
it: in *this* program it moves the 5-gon flip count by a factor of 2.3
(6,020 → 2,612) and by 3.0 on 4-gons (2,561 → 853), where closing the gap to the
original's figures would need a factor near 50. The remaining difference is some
further divergence in enumeration or filtering that prose cannot settle; putting
the two sources side by side is what would settle it, and **nothing in the ruling
depends on it**, because the quantity the ruling rests on is `SURVIVORS`, which
is zero in every block of both.

**The bound, stated honestly. This is evidence, not a theorem.** Two ring sizes,
three lattices, and the searched space is bounded by construction. It does *not*
cover n > 5, non-convex configurations beyond what 5 vertices reach, or rings
spanning many cells. One thing it covers more of than it appears to: the
parameter is the **step-to-spacing ratio**, not the spacing, because `snap` is
exactly scale-covariant on dyadic spacings — scaling the points and `h` together
changes only an exponent and leaves every grid index, every orientation sign and
every area sign identical. So "at `SnapGrid{1.0}`" is not the limitation it reads
as; "at 2, 4 and 8 sub-cell steps" is.

**Consequences, all three measured rather than argued.**

- The demotion of `RingDegenerateAfterSnap` to a self-check rests on measured
  flips with zero survivors, not on the structural sketch above.
- **Mutant 8 stays struck, for a measured reason.** There is no killing input in
  the searched space, and the rule stands: *a mutant with no killer is not a
  weaker mutant, it is a claim the round would have shipped as covered.*
- **`05-noder.md` risk 3's relocation is measured, not inferred.** Delete the
  winding re-check and every one of these rings is still refused, by
  `NonSimpleRing`, via 14(b) — which is exactly what the graze column counts.

**What would reopen it** is one survivor at any size, which promotes the status
back to a diagnosis and hands mutant 8 its killer. Re-running the program above
at a larger `n` or a wider lattice is the cheap way to look, and it is the right
thing to do before anyone *relies* on this beyond 5b. Note from the finer-step
block which direction to search: more vertices per cell is the wrong way, since
`repeat-id` eats the candidates before the question is asked.

**Reporting: one status, one message, and the message carries the count.** The
noder does not return a diagnostics vector, and the reason is
`include/terrain/cdt/result.hpp:12-17`'s distinction. `PslgBuilder` collects
because its input is authored outside the process and N broken chains are N
independent fixes. The noder's input is an *already valid* `Pslg`; every failure
above is snap-induced, and they share a single lever — the spacing. A vector here
would be a vector of the same fix repeated. What the vector would have carried
that a single status does not is the *extent*, so the message carries it:
`"3 rings collapsed to fewer than 3 nodes at spacing 0.5; first is chain 7"`.

## Degeneracy and failure policy

`05-noder.md`'s table is the input-condition table and is not restated. This is
the complementary one: **what each `NodeStatus` means, what raises it, and what
the person holding it should do.** The fourth column is the `CdtStatus` the same
input produced at increment 4, which is where the failure moved from.

| `NodeStatus` | Raised by | Actionable lever | Was, at increment 4 |
|---|---|---|---|
| `Ok` | — | — (but see "an open chain that collapses to one node", below) | `Ok`, or `NotNoded`, or `DegenerateGeometry` |
| `NotRun` | default-constructed `NodeOutcome` | none; a self-check | `CdtStatus::NotRun` |
| `InvalidSnapSpacing` | step 0, `is_valid_spacing` | supply a finite positive spacing | n/a |
| `CoordinateOutOfRange` | step 0, `can_snap`, naming the vertex | coarsen the spacing, or re-project | `Ok`, or nonsense |
| `RingCollapsed` | ring re-validation, naming the chain and spacing | **finer** spacing | `Ok`, or `DegenerateGeometry` |
| `RingDegenerateAfterSnap` | ring re-validation, `orientation<K>` | none; a self-check — see below | `Ok`, with a hole meshed as an island |
| `NonSimpleRing` | ring re-validation, repeated node id | polygon repair, upstream | `NotNoded` |
| `NotConverged` | `build<K>()` refused 14, at `max_rounds` | **finer** spacing, or raise the cap — **except on risk 16's corner graze, where neither helps** | n/a |
| `MalformedOutput` | `build<K>()` refused 11/12/13/16 | none; this is our bug | n/a |

**The lever for `NotConverged` has an exception and the message cannot tell you
which case you are in.** Risk 1 — snap rounding creating a crossing that was not
in the input — is genuinely cured by a finer spacing. Risk 16's corner graze is
not reliably cured by one: refining moves the lattice, so whether the flanking
cells still hold nodes is re-rolled rather than fixed, and a spacing that cures
one corner can create another elsewhere in the same data. For any *given*
fixture a different spacing usually does clear it, which is why the lever is
still worth printing; what it is not is a guarantee.
5b ships both under one status because the noder cannot distinguish them without
5d's predicate, and naming that here is cheaper than a status row that would have
to be deleted when 5d lands.

Note the shape `05-noder.md` names and it survives re-derivation: the failures do
not disappear, they **move**, and they move to the stage that can say something
actionable. Three of the nine rows point at the spacing, two of them in opposite
directions — `RingCollapsed` and `NotConverged` both want a finer grid,
`CoordinateOutOfRange` wants a coarser one. That is `05-noder.md` risk 2 showing
through the status enum, and it is why there is no default spacing to hide it.

### An open chain that collapses to one node is `Ok`, and it is pinned

At a coarse enough spacing every vertex of a short breakline snaps to the same
node. The chain survives as `count == 1` and contributes zero edges. The table
above had no row for it and the builder's guarantee list did not cover it, so it
was decided by whichever `if` happened to be written. It is ruled here.

> **Ruling: `node<K>` and `NodedPslgBuilder::build<K>()` accept an open chain of
> one node, with status `Ok`. The `count < 3` refusal is `is_closed(role)`'s
> alone, and that asymmetry is deliberate.**

Four reasons, in the order that decides it:

- **It is not our bug, and `MalformedOutput` is the only refusal available.**
  Every status the enum offers here reads as a diagnosis of the wrong thing.
  `MalformedOutput` says "this is our bug" — it is the self-check for guarantees
  the driver establishes *by construction* — and a chain collapsing under an
  over-coarse spacing is not a construction failure, it is the spacing doing
  exactly what the caller asked. `RingCollapsed` is the closed-chain form of the
  same event and its name would be a lie on an open chain. A tenth enumerator
  would cost a row, a `describe` arm, a Python enum member at 5c and a mapping,
  to report an outcome that is already correct.
- **Guarantee 16 forbids the alternative.** Dropping the chain would break chain
  order index-for-index with the input, which is the mitigation for half of
  `05-noder.md` risk 5 and the reason a Python-side per-chain attribute array
  stays valid without a map. A silently dropped chain is the failure 16 names.
- **The output is well-formed, not merely tolerated.** `edge_count(c)` is
  `count - 1 == 0`, `edge_base` still tiles, `edge_properties` is still
  index-aligned, and the node itself is in the vertex set and reaches the CDT as a
  point. Nothing downstream has to special-case it; the degenerate case is the
  empty range, which every loop already handles.
- **The closed case genuinely differs.** A ring of fewer than three nodes has no
  interior and no winding, so `orientation<K>` and the hole classification behind
  it have nothing to answer. An open chain of one node is asked no question it
  cannot answer.

**A caller who wanted the edge back wants a finer spacing**, which is
`RingCollapsed`'s lever arriving without `RingCollapsed`'s status. That the noder
does not say so is the accepted cost: it cannot distinguish "collapsed" from
"deliberately coarse" without a policy it has no business holding, and policy
lives at the composition root (`cli.py`, 5c).

**And it must be pinned, because an unpinned acceptance is a refusal waiting to
be tidied in.** The mutant that expresses "an open chain of one node ought to be
refused", spelled as minimally as a maintainer would spell it, is **`|| count < 2`
added to `noded_pslg_builder.hpp`'s guard**: it refuses `count == 1` and leaves
every other chain alone. That is the one that was unpinned. Two fixtures, both now
written, each the sole killer of its side:

- `tests/cpp/unit/test_noding_noded_pslg_builder.cpp` (`797144c`) — a
  **hand-built** candidate with an open chain of `count == 1`, required to
  `build<K>()` as `Ok`. This pins the guard's asymmetry at the guard, against
  arrays the noder did not produce, which is the whole reason the builder is a
  separate type.
- `tests/cpp/unit/test_noding_node.cpp` (`9600ad9`) — a short open breakline at a
  spacing coarse enough to collapse it, asserting `Ok`, `chains().size()` equal to
  the input's, and `edge_count(c) == 0` for that chain. This pins the end-to-end
  outcome and guarantee 16 together.

Measured here, not asserted, with the guard patched into
`include/terrain/noding/noded_pslg_builder.hpp:217`, the header `touch`ed after
each restore per `.claude/REQUIRED-READING.md`'s C++ stale-artifact hazard, and
`ctest --test-dir build` run over all 611:

| Mutant | ctest cases failed |
|---|---|
| baseline | 0 of 611 |
| `\|\| count < 2` added to the guard | **2** — exactly the two cases above |
| `is_closed(ch.role) && count < 3` widened to bare `count < 3` | **26**, across four suites |

**An earlier revision of this paragraph named the wrong mutant and its premise
was false.** It said "nothing currently fails if someone widens ... to
`count < 3`". That widening is killed 26 times over, because every open chain in
the builder suite's `valid_candidate()` has `count == 2` and the widened guard
refuses them all — it was never the unpinned one. The figure above is the one
measured for this section; `@tester` reported **8** for the same mutant and
`@orchestrator` **26**, and the units differ (assertions against registered
`ctest` cases) rather than the measurements disagreeing. **No reconciliation is
attempted, deliberately: the claim was "nothing fails", and that is refuted at
8, at 26, and at 1.** Recording which figure is which unit is worth more than
picking one.

**The conclusion is unchanged.** The acceptance genuinely was unpinned — at the
commit before those two fixtures, the `count < 2` mutant killed nothing, which is
visible in the row above: the only two cases it fails are the two that section
added.

**And the closed arm of the same guard had no test at all.** An earlier revision
of the first bullet called the closed-chain-under-3 refusal "the same file's
existing" case and asked that the two read as a pair. There was no such case:
`is_closed(role) && count < 3` was reachable, argued for in this document, and
unexercised, so deleting the guard outright was as free as widening it. `@tester`
wrote it in `797144c` alongside the acceptance case, and it is the sole killer of
guard deletion. **Nothing in the guarantee list above ever covered it** — 11-16
say what a `NodedPslg` promises, not which refusals have fixtures, and reading a
guarantee as a claim about coverage is how this gap stayed open. The pair exists
now; what did not exist is the half this document assumed.

### What becomes unreachable in `CdtStatus`, at 5c

Unchanged from `05-noder.md`'s ruling, recorded here with the file that has to
change:

- `PointOnConstrainedEdge` and `ConstrainedEdgeIntersection` — guarantee 14 makes
  both unreachable, so `CdtStatus::NotNoded` joins `MalformedInput` as a
  self-check.
- `DuplicatePointsFound` — guarantee 12 plus `world` injectivity.
  `PolylineDuplicateConsecutivePoints` — guarantee 13. `DegenerateGeometry`
  loses both its reachable rows; `AllPointsAreCollinear` was already unreachable.
- **Reachable from a `NodedPslg`: `Ok`, `InvalidTopology`, `BackendFailure`.**

**This is 5c's documentation change, not 5b's**, because it is only true once the
signature has changed. It lands as one paragraph appended beneath `04-cdt.md`'s
mapping table, leaving the table itself intact as the record of what was true at
increment 4. `describe(CdtStatus::NotNoded)`'s text
(`include/terrain/cdt/result.hpp:46-48`) keeps naming the noder as the fix,
because at that point it is naming the thing that should have been called.

### Nesting is still not computed

`03-pslg.md:520-526` reserved `nesting_forest<K>(const NodedPslg&)` for the
increment where non-crossing input makes the answer meaningful, which is this
one. Still deferred, and for `05-noder.md`'s better reason rather than the
original one: `detria` computes nesting itself and diagnoses the failures as
`InvalidTopology`, so a second implementation here is derived data with no unique
consumer. Named a third time so it is not re-derived a fourth.

## What the signature change actually costs

5c, in full, because `05-noder.md` books it as one ~10-line row.

**C++, and this part really is mechanical.** `include/terrain/cdt/triangulate.hpp`
(the `CdtBackend` concept and the generic entry point),
`include/terrain/cdt/detria_backend.hpp` (the declaration) and
`src/cdt/detria_backend.cpp` (the definition) each replace `const Pslg&` with
`const NodedPslg&` and swap the include. The wrapper body is untouched: it walks
`vertices()`, `chains()` and `indices_of(c)`, all of which `NodedPslg` has with
the same signatures. ~6 lines.

**`bindings/core.cpp`, and this part is not.** `bindings/core.cpp:466` calls
`terrain::cdt::triangulate<DetriaBackend>(pslg, options)` with a `const Pslg&`;
retyping the C++ entry point breaks that translation unit, so the binding is
**not optional and cannot be deferred**. There is no smaller consistent step
than: `py::enum_<NodeStatus>` plus a `describe` overload, `py::class_<SnapGrid>`,
`py::class_<NodedPslg>` mirroring the five `Pslg` properties plus `grid`,
`edge_properties` and `node_of_input_vertex`, `py::class_<NodeOutcome>`, and an
`m.def("node", ...)` releasing the GIL the way `triangulate` does. Python cannot
otherwise obtain the type the CDT now demands.

**`src_python/tin_engine/_core.pyi`** gains stubs for all of it or `mypy --strict`
fails on `cli.py`'s first call.

**`src_python/tin_engine/cli.py`.** `_triangulated` gains a `node()` call between
`build_pslg` and `triangulate`, a `--snap-spacing` option with a default *in
Python* (the composition root is where policy lives; the C++ ruling is
unaffected), and a second failure presentation — a noder failure has a
`NodeStatus` and a message, not a `CdtStatus`, and the header band must carry the
noder's own words the way it already carries the backend's.

**And the scene must be fed the `NodedPslg`**, per "What the renderer buys this
increment". That is a one-line change in `cli.py` and it is the one line that,
got wrong, turns the whole gallery red while the mesh underneath is correct.

**Python suites change too** — `tests/python/test_core_cdt.py`'s
`test_reports_a_non_noded_input_as_a_failure_status` asserts that
`crossing_pslg` fails, and after 5c that input cannot reach `triangulate` at all.
No line number is cited, because 5c is the increment that moves it; locate it
with `grep -rn non_noded tests/python/`, and see `05c-noder-wiring.md`'s churn
list for what it became — it split rather than moved, and the CDT half is a
`TypeError` rather than a status.
That is `@tester`'s change, in its own commit with the reason in the message, per
`docs/increments/README.md`'s "no test change hides inside an implementation
commit". Tests are excluded from the LOC count; the *round* is not free.

## Files and LOC

**Two instruments, declared, because the repo has two and they disagree.** For
C++, `grep -vcE '^\s*(//|$)'` — preprocessor directives are code. For Python,
`grep -vcE '^\s*(#|$)'`. `06-cdt-viewer.md`'s single combined command undercounts
C++ by every `#include`; see premise 3. Neither instrument excludes docstrings or
pybind `R"doc(...)"` bodies, and `CLAUDE.md` §2's unit — non-comment lines — is
the definition; these are the instruments that approximate it.

Re-runnable for 5a, the calibration point:

```sh
for f in include/terrain/core/snap_grid.hpp include/terrain/noding/intersect.hpp \
         include/terrain/noding/node_set.hpp; do grep -vcE '^\s*(//|$)' $f; done
```

prints `57`, `89`, `33` — 179, which is the figure `41ac5b6` reports.

**5b — the noder in C++:**

| File | Contents | Est. LOC |
|---|---|---|
| `include/terrain/core/noded_pslg.hpp` | `NodedPslg`, its accessors, `edge_base` | ~85 |
| `include/terrain/noding/broad_phase.hpp` | bbox-bucketed segment index, `for_each_candidate(Box2, F&&)` | ~70 |
| `include/terrain/noding/noded_pslg_builder.hpp` | `NodeStatus`, `describe`, `NodeOutcome`, `NodedPslgBuilder::build<K>` | ~110 |
| `include/terrain/noding/node.hpp` | `NodeOptions`, `node<K>`, the six-step pass, ring re-validation | ~180 |

**~445 production LOC.** Header-only; nothing in `src/noding/`, because the
driver is a template on `K` exactly as increment 3's validator is. Non-C++
changes in the same PR, listed separately because they are not production code:
**five** suite registrations in `tests/cpp/CMakeLists.txt` — four headers but
five files, because `node.hpp` has the second file declared under "What is worth
testing" — and they are **not all registered the same way**.

Four go through `add_terrain_backend_test`, because they name a kernel.
`prop_noding_broad_phase` goes through plain `add_terrain_test`
(`tests/cpp/CMakeLists.txt:207`), **and the asymmetry is a design property, not
an oversight.** `broad_phase.hpp` includes `core/bbox.hpp`, `core/point.hpp` and
`core/segment.hpp` and no predicate header: it compares bounding boxes and never
asks an orientation. Registering its suite without the backend leaves
`terrain_predicates` off that link line, so a broad phase that ever reaches for a
predicate **fails to link** rather than passing quietly. Under
`add_terrain_backend_test` the same regression would compile, link, and leave
every assertion green — the proof deleted and nothing red. The kernel-free
one-directional bbox contract at the head of the header is a claim about the
*link line* as much as about the code, and this is the registration that checks
it. Anyone tidying the four-versus-one into uniformity is removing the check.

**5c — the crossing:**

| File | Contents | Est. LOC |
|---|---|---|
| `include/terrain/cdt/*`, `src/cdt/detria_backend.cpp` | `const Pslg&` → `const NodedPslg&` | ~6 |
| `bindings/core.cpp` | `NodeStatus`, `SnapGrid`, `NodedPslg`, `NodeOutcome`, `node()`, retyped `triangulate` | ~215 |
| `src_python/tin_engine/_core.pyi` | stubs for all of the above | ~90 |
| `src_python/tin_engine/cli.py` | the `node()` call, `--snap-spacing`, the noder's failure presentation | ~30 |

**~341 production LOC.**

**Together ~786, which is why they are two PRs.** The binding row is estimated at
~215 rather than the ~108 a line-by-line sketch suggests, because every
binding estimate this project has made has come in at roughly twice: 6a's
`bindings/core.cpp` at +283 against ~150 and `_core.pyi` at +125 against ~60.
Applying the measured factor rather than the sketch is the whole of Gate C's
reason for existing. The C++ rows are *not* doubled, for the complementary
measured reason: increments 3, 4 and 5a all came in at or under their C++
estimates, 5a by a third.

**The edge-property correction does not move either figure**, and that is worth
recording rather than asserting. Re-checked row by row: `noded_pslg.hpp` swaps
one span's element type; `node.hpp`'s merge stays one `|` per contributor and
one array; `noded_pslg_builder.hpp` does not check guarantee 15 at all;
`bindings/core.cpp` exposes a `py::array_t<std::uint32_t>` where it would have
exposed a `py::array_t<std::uint8_t>`, at identical length. **5b stays ~445 and
5c stays ~341.** The type the noder merges is not the noder's to define, which is
exactly why widening it costs the noder nothing.

**Increment 7, sized here only so nobody assumes it is small or assumes it is
large**: `core/edge_properties.hpp` ~45, `core/pslg.hpp` and
`core/pslg_builder.hpp` ~10 between them, `bindings/core.cpp` ~25,
`_core.pyi` ~12, `viz/protocols.py` ~4, `viz/scene.py` ~10, and on the
`increment6b-ii-renderer-red` branch `viz/style.py`, `viz/svg.py` and
`viz/fixtures.py` ~25 — **~130 production lines**, comfortably one PR. The round
is not cheap for a different reason: the test churn crosses eight files and two
of them are invariant-critical, so the mutation budget is spent on suites that
already have one.

**The estimate survives the red round, and one row is the one to watch.**
`@tester` built a throwaway implementation of all four headers during the red
round and measured it, under this section's C++ instrument, at **350** non-comment
lines — 95 under the estimate, with only `noded_pslg_builder.hpp` over its own
row, at 127 against ~110. That reading is confirmed: a throwaway is a *lower*
bound on the shipped header, since it carries no comments the project's density
would demand and none of the header's documentation burden, so 350 against 445 is
headroom rather than a refutation, and the one row that came in over is the row
whose contents this round has grown. None of this round's rulings adds production
lines — 5d is deferred, the `NodeStatus` reading is a comment, and
`RingDegenerateAfterSnap`'s demotion deletes no check — so **5b stays ~445**.

### Reconciliation: ~445 estimated, 675 measured

Measured at **the last commit on this branch that touches a production file**.
The anchor is a rule and not a hash, for `07-edge-properties.md`'s reason: a hash
names a commit, the commit that updates the hash is itself a commit, and a
correction round that touches a header falsifies the line it just wrote. Resolve
it:

```bash
git log --oneline -1 -- include/ bindings/ src_python/
```

and re-run this section's C++ instrument per file:

```sh
for f in include/terrain/core/noded_pslg.hpp include/terrain/noding/broad_phase.hpp \
         include/terrain/noding/noded_pslg_builder.hpp include/terrain/noding/node.hpp; do
  grep -vcE '^\s*(//|$)' $f; done
```

| File | Est. | Actual | |
|---|---|---|---|
| `include/terrain/core/noded_pslg.hpp` | ~85 | 72 | under |
| `include/terrain/noding/broad_phase.hpp` | ~70 | 84 | |
| `include/terrain/noding/noded_pslg_builder.hpp` | ~110 | 249 | **2.3×** |
| `include/terrain/noding/node.hpp` | ~180 | 270 | 1.5× |
| **Total** | **~445** | **675** | **1.52×** |

**The row this section named as the one to watch is the row that blew out**, and
that is the finding rather than the total. The paragraph above says "the one row
that came in over is the row whose contents this round has grown", naming
`noded_pslg_builder.hpp` at 127 against ~110 in the throwaway. Shipped, it is 249
— 139 over its own row, and 60 % of the whole overrun. The warning was correct in
direction and understated by an order of magnitude in size, which is the more
useful record: a row identified as growing was still estimated as if it were not.

**Where the 230 went, per row rather than in aggregate**, because
`07-edge-properties.md`'s own correction round was forced by a total nobody had
decomposed. `noded_pslg_builder.hpp` carries six guarantee checks as six named
private members, each with its refusal message formatted and its ordering reason
in a comment, plus `NodeStatus`, `describe` and `NodeOutcome`; the estimate priced
the checks and not the diagnostics — a `std::format` refusal that names the chain
and the spacing is three or four lines where "return an error" is one, and the
status table above owes one such message per row. `node.hpp` at 270 against ~180
is the six-step pass plus `revalidate_rings` as a separate function, which the
estimate folded into the step list. `broad_phase.hpp`'s +14 is the documented
one-directional contract; `noded_pslg.hpp` came in under.

**The instrument counts the comments this project asks for**, and that is not an
excuse but it is an accounting fact worth carrying: `grep -vcE '^\s*(//|$)'`
excludes `//` and blank lines, so the non-comment count is what `CLAUDE.md` §2
means and these headers' unusually heavy `//` blocks are already excluded. The
overrun is code, not prose — the opposite of `features.py`'s at increment 7, and
the reason to state which of the two a given overrun is.

**The ceiling holds and is not close**: 675 against 700 non-comment production
lines, 25 of headroom. That is thin enough to be worth saying out loud — **5b has
no room for another production line**, and any ruling that would add one belongs
in 5c or 5d rather than here. Every ruling this round added is prose or a test.

**What this does to 5c's ~341.** The C++ estimate's track record is now
**three-for-four, not three-for-three**: increments 3, 4 and 5a came in at or
under, and 5b came in at 1.52×. The sentence above that declines to double the C++
rows — "for the complementary measured reason" — rests on that record, so it is
weaker than when it was written. It is not withdrawn, because 5c's C++ row is ~6
mechanical lines and the estimate's risk is concentrated in the binding row, which
is already doubled on the measured binding factor. What changes is the margin:
~341 was quoted against 700 as comfortable, and the honest reading now is that a
5b-sized 1.5× on the Python and binding rows would put 5c around 500 — still one
PR, no longer comfortable, and inside one row's slip of **Gate C's 550**. Gate C
already obliges 5c's design to re-estimate `bindings/core.cpp` and `_core.pyi`
against what 5b's types turned out to be; this table is the evidence that the
obligation is real, and it now applies to the C++ rows too.

**What to carry into the next estimate.** The row that grew during the design
round is the row to re-estimate, not merely to watch: this document flagged
`noded_pslg_builder.hpp` as having grown and then left its ~110 unchanged. Naming
a row as at risk and not re-pricing it is the same move as an estimate nobody
reconciles.

**Nobody should relitigate the seam on a line count in either direction.** If 5b
comes in at 320 the seam still stands, because 5c's Python surface would not fit
beside it, and because the two halves are reviewed by different eyes against
different gates — a C++ suite under ctest and sanitizers, a Python suite under
mypy, ruff and an 85 % coverage floor.

## What is worth testing

**Four suites, one per production header — Gate B's mapping.** Three carry a
mutation round; one deliberately does not, and the reason is given rather than
left as an omission (`05-noder.md`'s precedent: `node_set.hpp`'s suite was
excluded in writing).

| Header | Suite | Invariant-critical? |
|---|---|---|
| `broad_phase.hpp` | `tests/cpp/property/prop_noding_broad_phase.cpp` | **yes**, mutation round |
| `noded_pslg_builder.hpp` | `tests/cpp/unit/test_noding_noded_pslg_builder.cpp` | **yes**, mutation round |
| `node.hpp` | `tests/cpp/property/prop_noding_no_crossings.cpp` | **yes**, mutation round |
| `noded_pslg.hpp` | `tests/cpp/unit/test_noding_noded_pslg.cpp` | **no** |

`tests/cpp/unit/test_noding_node.cpp` is where the driver's hand-built fixtures
live — the crossing, the T-junction, the lattice-corner tie, the road-along-a-
river, each failure status — and it is `node.hpp`'s *second* file. That is the
one place Gate B's one-file-per-header rule is relaxed, deliberately and in
writing: a property file and a fixture file for the same header are a split
along a different axis than the one Gate B protects, and both can move to a 5b-ii
together if a seam is ever cut.

**What the lattice-corner fixture is for, stated explicitly because the finding
changed it.** It was written as a fixture that would confirm guarantee 14; it is
now the fixture that refutes it, so its purpose has to be in the design rather
than inferred from the assertions. It carries **two cases and they are
independent by construction**:

1. **The arithmetic, in 5a's predicates alone.** `segment_meets_cell<K>` on the
   diagonal against `(1,0)` and `(0,1)`, and on the anti-diagonal against `(0,0)`
   and `(1,1)` — the probe above, as a test. It calls nothing in 5b. **It stands
   whatever `node<K>` does**, including after 5d changes the predicate, at which
   point **all four** of its `segment_meets_cell` expectations flip to `false` —
   both pairs are grazes through the *same* lattice corner, so a predicate change
   that touches one touches all of them — and it becomes 5d's regression test
   without being rewritten. An earlier revision said "the second pair"; the count
   was wrong and the load-bearing half, *without being rewritten*, was not. That independence is the point:
   it pins the geometric fact that every other ruling here rests on, at a layer
   that no 5b decision can move.
2. **The outcome, `NotConverged`.** The full driver on the tie, asserting the
   status and that `outcome.pslg` is disengaged. This one is 5d's to change.

**And it is load-bearing either way**, which is the argument for carrying both
rather than folding them. Case 2 bears on mutant 5 — the loop condition spelled
"until a round produces no new node". Under that spelling the tie converges to
`Ok` at the end of round 1, because `(1,0)` and `(0,1)` are already nodes — they
came from the other chains — so round 1 splits an edge at two **existing** nodes
and manufactures none; under the correct spelling it runs to the cap.

**The precise claim, because a looser one contradicts the mutants list.** Mutant
5 below names the road-along-a-river fixture as its killer, and that fixture
exists; so case 2 is *not* "mutant 5's only killer", and an earlier revision of
this paragraph said it was. What is true, and what the suite records, is
narrower: **case 2 is the only fixture in the increment on which the two loop
spellings produce different *statuses* rather than different edge counts** — so
it is the only killer that fails loudly, and the others fail by a count a reader
has to interpret. That is the reason not to delete it as "a test of a bug we are
going to fix", and it is stated narrowly on purpose: a warning that rests on a
claim this document contradicts three screens further down is a warning someone
will argue their way past.

**`prop_noding_broad_phase.cpp` — invariant-critical, mutation round.** This is
`05-noder.md` risk 12's whole mitigation and the one function in the increment
with an oracle that does not touch the noder. For generated segment sets and
generated query boxes, the visited set is compared against **all-pairs bbox
intersection computed without the index**. The contract mentions no segment pair,
no node and no snap grid, which is exactly what makes that oracle possible; a
check written in terms of "the pairs the noder needed" would share the defect.

**`test_noding_noded_pslg_builder.cpp` — invariant-critical, mutation round.**
The only suite that can fail when the *verifier* is broken. It hands
`NodedPslgBuilder` candidates constructed **by hand, not by the noder**, and
requires a refusal for each: a live T-junction (14(b)), two crossing edges
(14(a)), a partial collinear overlap (14(a), the amended clause), two equal
vertices (12), an off-grid coordinate (11), a zero-length edge (13), a chain
slice that overlaps its neighbour (16). And, symmetrically and just as
important, an **acceptance**: the road-along-a-river duplicate-edge candidate,
which the pre-amendment 14(a) would have refused.

**A partial collinear overlap cannot be isolated to 14(a), and the suite must
not pin the clause.** `@tester` measured it: the overlap's inner endpoint is a
node that lies *on* the other edge, so its cell meets an edge it is not an
endpoint of, and 14(b) refuses it too. It is always a double refusal. The case
is still worth carrying — it is the amended clause's negative half — but the
assertion is on the **status**, which is `NotConverged` for either clause, and
never on which clause fired. Nothing in the design lets a caller see the clause,
so a suite that could see it would be reading the implementation. Same for the
crossing-edges case, which 14(b) also catches once the crossing point is a node.

**`prop_noding_no_crossings.cpp` — invariant-critical, mutation round.** The name
`testing.md`'s layout block reserves; `grep -n prop_noding_no_crossings testing.md`
resolves it, and the line number is deliberately not written down here, per
`.claude/REQUIRED-READING.md` — this citation was already stale once, having been
written as a bare line number that the layout block had since moved off.
It runs the full noder on generated constraint sets of
a few dozen segments, then checks, **with no broad phase anywhere in the check**:

- **Guarantee 14, brute force** over every (edge, edge) and every (node, edge)
  pair. A broad-phase defect and a matching verification blind spot cannot cancel
  against an oracle that does not share the generator. O(n²) is why it is
  small-input and a property rather than the production path.
- **Guarantee 15, from the input**, in the relation the producer used — the
  `segment_meets_cell` + arc-order-betweenness form `05-noder.md` fixes, never
  the dedup's provenance map, and never exact incidence. Stated over sets: for
  every output edge `e`, the union of the property sets of every input chain that
  contributes to `e` under that relation is a **subset** of
  `edge_properties()[e]`. **One direction only — subset, not equality**: the
  converse does not hold and must not be asserted, because a chain running
  *near* both nodes and spanning them satisfies the relation without having
  contributed, so the output legitimately may be a strict superset of what the
  oracle can prove. A spuriously *added* property is therefore not caught here
  and is owed a mutant instead — mutant 7 below.

  Note what the widening buys the oracle for free: under one bit, "the set is a
  superset" and "the bit is set" were the same assertion, so a merge that
  **replaced** rather than unioned was indistinguishable from a correct one
  whenever the replacing chain happened to be the river. Over sets the two come
  apart, and mutant 7 is killed by an assertion rather than by a fixture chosen
  to make the difference visible.
- **Guarantee 16**, chain-for-chain against the input.
- **Determinism**: the same input shuffled in vertex order produces node ids that
  are a function of the point set, hence identical output geometry.

`tests/cpp/property/noding_generators.h` arrives here — `testing.md` describes it
as producing "random sets of polylines with controllable density of
intersections", which is 5b's generator and not 5a's.
`grep -n 'controllable density' testing.md` resolves the line; it is not written
down, for the same reason as the citation two paragraphs above, and that one was
stale too when this round found it.

**`test_noding_noded_pslg.cpp` — not invariant-critical, and no mutation round.**
It is an accessor suite over a value type: spans, slices, `edge_base`'s prefix
sum, `ring(c)` not throwing, `const` concurrent read. Its failure mode is a typo
and it is loud. The README's rule is to spend the mutation budget where the
topology decisions are; none of them are here. The one accessor with arithmetic
in it, `edge_base`, is checked against `edge_count`'s running sum in the property
suite anyway, which is a stronger check than a mutant.

### Mutants the round must kill

Continuing `05-noder.md`'s numbering conceptually but restarted, since these are
a different increment's:

1. **`for_each_candidate` bucketing by segment endpoints** instead of by every
   bucket the bbox overlaps. Killed by a generated segment longer than a bucket
   whose interior crosses a query box its endpoints miss. This is risk 12's live
   defect and the reason the broad-phase property exists.
2. **A `Point2` overload of `for_each_candidate`, or the node query built from
   `world(g)` instead of `cell_min`/`cell_max`.** Killed by a node whose cell
   straddles a bucket boundary. The overload must not exist; if it appears in a
   diff, that is the mutant as a reviewable change.
3. **The split pass driven by `classify<K>` instead of `segment_meets_cell<K>`** —
   the provisional design 5a's audit refuted. Killed by the T-junction fixture,
   which classifies `Disjoint` at a decimal spacing and so is never split.
4. **`segment_meets_cell` queried only for the nodes of the segment's own
   chain.** Killed by the crossing fixture: the other chain's nodes are the ones
   that matter.
5. **The loop condition spelled "until a round produces no new node"** instead of
   "until the verification passes". Killed by a fixture whose second round splits
   an existing edge at an *existing* node — the road-along-a-river at a spacing
   where the road's vertices snap into the river's cells on round 1 and the
   river's edges must be split on round 2. **The lattice-corner fixture's case 2
   kills it too, and differently**: there the two spellings differ in *status*
   rather than in edge count, which is the only loud failure of the two. See
   "What the lattice-corner fixture is for" — the two entries are consistent and
   were not always.
6. **`max_rounds` treated as a success** — returning the last candidate instead of
   `NotConverged`. Killed by asserting that a non-`Ok` outcome's `pslg` is
   disengaged, on a fixture built to need more rounds than a cap of 1.
7. **The property union replaced by "take the first contributor"** or by "take
   the chain's own set". Killed by the road-along-a-river fixture in
   `test_noding_node.cpp`, asserting the merged entry on the *road's* flat
   position and requiring it to contain **both** properties — which is the
   assertion the one-bit form could not make, since under it the road's bit was
   the absence of the river's. This is the mutant `05-noder.md` says 5b owes
   because guarantee 15's oracle is one-directional, and it is the only thing
   standing between a lost property and a silent wrong answer.
8. **Struck, and the reason is the finding rather than a change of mind.** It
   read: "the winding re-check dropped, or a reversed ring silently repaired;
   killed by the sliver-hole fixture whose winding flips under snapping". There
   is no such fixture — see "`RingDegenerateAfterSnap` is a self-check". Four
   candidates were measured and all four are refused as `NonSimpleRing` by the
   split pass before the winding check runs, so this mutant has **no killer
   input**, which is exactly what it means for the check to be a self-check.
   Its live half is mutant 10, which kills the same silent wrong mesh through
   the node repeat that actually occurs. **A mutant with no killer is not a
   weaker mutant; it is a claim that the round would have shipped as covered.**
   The bounded ring search in that section has since measured it: across two
   ring sizes, three lattices and two independently written programs, **every
   winding flip grazes and no block has a survivor**. So the absence of a killer
   is measured rather than assumed. Reinstated the moment that search turns up a
   survivor at any size.
9. **The `RingCollapsed` check dropped.** Killed by a hole smaller than one cell.
10. **`NonSimpleRing` checked by a segment-pair scan that only finds interior
    crossings**, missing a ring that revisits a node. Killed by a figure-eight
    ring whose crossing is *at* a node after noding — which, by 14(a), is the only
    form it can take.
11. **Arc order by node id instead of by dot product.** Killed by a segment whose
    split nodes sort differently under `GridPoint`'s lexicographic order than
    along the segment. An earlier revision said "any segment running
    down-and-right"; that is wrong and would have produced a mutant that
    survives. `GridPoint`'s order is lexicographic in `(ix, iy)`
    (`include/terrain/core/snap_grid.hpp:43-51`), so the direction that inverts
    it against arc order is **decreasing `ix`** — down-and-right has `ix`
    increasing and the two orders agree. The fixture is a segment whose `ix`
    decreases from `a` to `b`, or a vertical segment whose `iy` decreases.
12. **An input-dependent anchor**: subtracting a bounding-box origin before
    snapping. `05-noder.md` mutant 12 assigns this to 5b's driver, since 5a's
    `snap` has no argument it could arrive through. Killed by noding a fixture,
    then noding it again with one extra far-away unreferenced vertex.

    **The assertion is not "bit-identical output", and getting that wrong costs
    the mutant.** Step 3 builds the `NodeSet` over **every snapped input vertex**,
    referenced or not — deliberately, because that is what makes
    `node_of_input_vertex` total over the input's vertex array — so the extra
    vertex *is* a node and the second run has exactly one more. Node ids are the
    set's sorted order, so the extra vertex also shifts every id that sorts after
    it, and a naive index-by-index comparison fails on correct output. Two
    assertions restore it, and the second is a property of the fixture rather
    than of the code: the node count is exactly one greater, and the extra
    vertex is **placed so that it sorts last** — a grid index larger in `ix`
    than anything in the domain — so no existing id moves and the
    index-by-index coordinate comparison is legitimate again. Placing it there
    is part of the fixture, not an accident of it. The mutant dies because an
    anchored `snap` moves the bounding box and therefore moves *every* node, so
    the comparison fails at the first index rather than only at the last.
13. **The verifier's clause 14(b) returning `true` unconditionally**, or checking
    `on_segment<K>` instead of `segment_meets_cell<K>`. Killed **only** by
    `test_noding_noded_pslg_builder.cpp`'s hand-built T-junction candidate — not
    by any property over the noder's own output, which is the argument for the
    builder being a separate type.
14. **The `can_snap` admission pass dropped.** Killed by a coordinate past
    `kMaxGridIndex`, which must be `CoordinateOutOfRange` and not an assertion
    failure or, in release, a wrapped index.

**Not mutants, listed so nobody spends the round on them.** The edge-key dedup
keyed on `world()` coordinates rather than on node ids: behaviourally invisible
under `world`'s injectivity, exactly as `05-noder.md`'s mutant 11. And the broad
phase's bucket *count*: any positive count is correct, only slow, so no fixture
can distinguish `sqrt(n)` from 1 except by timing, and a timing assertion is not
a test.

### Template spend

**Zero `TEMPLATE_TEST_CASE`. Everything runs under `DefaultKernel`.**

`broad_phase.hpp` is kernel-free — it compares bounding boxes and never asks an
orientation. Everything kernel-dependent in `node.hpp` and
`noded_pslg_builder.hpp` reaches `K` through `classify<K>`,
`segment_meets_cell<K>` and `orientation<K>`, and that those flow through `K`
rather than through inline doubles is what 5a's single cross product already
proves (`05-noder.md`, "Template spend": the near-parallel crossing at Web
Mercator magnitudes). A second cross product here would re-prove 5a's point at
the cost of compiling the two heaviest files in the increment twice. The README's
rule is that a cross product must prove something a later increment depends on;
this one would not.

**`FastKernel` is forbidden in the property suites**, same ruling and same reason
as increments 4 and 5a: the oracle is built from the kernel, and an oracle less
exact than its subject reports false failures on precisely the degenerate inputs
the property is about. `05-noder.md`'s second reason applies with more force
here than anywhere — snapped data is where the filter fall-throughs live, 38 % of
grid-collinear triples at 0.1 m, and every input in this increment is snapped by
construction.

## Risks

Continuing `05-noder.md`'s register, which ends at 12.

13. **The signature change is not mechanical, and three documents say it is.**
    `03-pslg.md:490-493` and `04-cdt.md:659-660` both book it as mechanical and
    `05-noder.md` estimates it at ~10 lines; `bindings/core.cpp:466` is the line
    that refutes all three. Mitigated by making it 5c's own PR with the Python
    surface it drags along. Residual: a reader of any of those three documents
    still arrives expecting a one-line change, which is why the correction is
    made in `05-noder.md` itself in this PR rather than recorded here.
14. **The binding estimate is the one most likely to be wrong, and it has
    doubled twice.** 6a: `bindings/core.cpp` +283 against ~150, `_core.pyi` +125
    against ~60. Mitigated by estimating at the measured factor and by Gate C.
    Residual: if it doubles *again* from the already-doubled figure, 5c splits.
15. **The scene join is index-based, so 5c has exactly one line that must feed
    the `NodedPslg` and not the `Pslg`.** Getting it wrong produces a picture
    that is alarming and a mesh that is correct — the most expensive kind of
    wrong, because the instrument is the thing that lies. Mitigated by
    `tests/python/test_viz_protocols.py`'s join, which is where `PslgLike`
    conformance is checked, and by the `not-noded` fixture rendering clean.
16. **An exact corner graze makes the split pass and guarantee 14(b)
    inconsistent, and the input cannot converge.** Upgraded from "named" to a
    ruling, because the worked case above was checked and came back red; the
    probe and the divergent orbit are in "The tie does not converge" above and
    are not repeated here.

    **What actually triggers it, corrected.** An earlier revision of this entry
    said "on axis-aligned cadastral data at a decimetre spacing that is not
    exotic". The probe's lines 3 and 4 refute that: an axis-parallel edge grazes
    no flanking cell, and axis-aligned runs are the *clean* case. The trigger is
    a constraint edge passing **exactly through a lattice corner** — a 45° run,
    or any slope whose grid crossings land on corners — with nodes present on
    **both** flanking cells. A 45° parcel or property boundary at a decimetre
    spacing is not exotic either, so the frequency claim survives its own
    correction; only the geometry in it was wrong. It also does **not** require a
    constructed point: `docs/increments/05d-corner-graze.md` carries a
    configuration of round decimetres on a decimetre grid where `classify<K>`
    reports `Disjoint` for every pair and the defect appears anyway.

    **Ruling: `NotConverged` is the right answer for 5b and the wrong answer for
    the product, and the fix is `docs/increments/05d-corner-graze.md`**, owned by
    `@architect`. Not 5b: the fix changes 5a's shipped `segment_meets_cell`,
    which has its own mutation round and mutants 13–17 aimed at it, and
    reopening it inside a PR that is already at ~445 lines and whose suites are
    committed red buys a second unreviewed decision for the price of one. Not 5c
    either: 5c is the crossing, its budget is Python surface, and a predicate
    change landing in the same PR as the binding would be invisible under it.

    **The specification is in that file and is not duplicated here** — what the
    replacement predicate is, why it is not half-open cells, which existing
    assertions invert, and what it costs the 5a mutation round. This entry states
    only what 5b ships in the meantime: **`NotConverged`**, with the message
    naming the spacing, plus the lattice-corner fixture, which is what makes the
    defect visible the moment anyone re-opens the question. The residual is that
    between 5b and 5d a user with a 45° breakline can be told to refine a spacing
    that may not help them. Refining does not remove the pathology, it re-rolls
    it: the lattice moves, and whether the flanking cells of some corner on that
    run still hold nodes is a fresh question at every `h`. That is the part a
    reader must not miss, and it is why the lever column for `NotConverged` in
    the status table carries a second sentence — which 5d deletes.

17. **The verifier shares `segment_meets_cell<K>` with the driver.** That is
    borrowing a *predicate*, which `computational-geometry/SKILL.md` licenses,
    and not borrowing *records*, which it forbids — but it is one shared object,
    and if `segment_meets_cell` were wrong the split and the verification would
    be wrong together. Mitigated one layer down: `segment_meets_cell` is 5a's,
    with its own mutation round and five mutants (13–17) aimed at exactly this.
    Residual: 5b inherits the quality of that round and cannot improve on it.
18. **5b ships with no production caller, for the second increment running.**
    Only tests instantiate `node<K>` until 5c. `05-noder.md` names this as 5a's
    honest cost; it is now a two-increment cost, and the sequencing risk
    (increment 3 risk 1, increment 4 risk 4, `05-noder.md` risk 6) persists
    through both. Discharged at 5c and not before.
19. **`prop_noding_no_crossings` is O(n²) per round and CI runs it twice** —
    ubuntu and macos, plus an asan+ubsan Debug build of the same suites. At a few
    dozen segments and a handful of rounds that is milliseconds; the risk is that
    "a few dozen" drifts upward during the round because a bigger generator finds
    more. Mitigated by fixing the generator's size in the suite with a comment
    saying why it is small. Residual: a defect that only appears above that size
    is not found here. The gallery is the counterweight, and it is the third
    thing the renderer buys.
20. **The edge-property vocabulary is Python's, and no C++ suite can check
    it.** The core merges opaque bits; the bit-to-name mapping lives at the
    boundary, so two producers that disagree about which bit means "river"
    produce a mesh that is wrong in a way nothing in `terrain::` can see.
    Mitigated at the boundary and not here: the mapping is one Pydantic model
    rather than a convention, with a digest over its `(bit, name)` pairs for an
    artifact to carry and compare on read
    (`src_python/tin_engine/features.py`; `docs/increments/07-edge-properties.md`
    inherits this entry as its own risk and owns the mitigation). **5b still
    cannot check it**, because the noder sees bits and never reads a
    fingerprint. Residual and unchanged: this is the price of keeping feature
    names out of `terrain::`, and it is the right price, but it is a price.
21. **`max_rounds` has a default and a default is a policy.** Four is defended
    above on measured rarity, not on proof. If real breakline data needs more,
    the symptom is `NotConverged` on input a person believes is fine, and the
    lever is a `NodeOptions` field rather than a rebuild — which is why it is a
    field.

## Documentation this PR fixes

Per `docs/increments/README.md`: fixed here, or not recorded. There is no ledger.

**Six of `05-noder.md`'s eleven promised fixes did not land in the 5a PR**, and
that is itself the defect this section opens with. Re-runnable:

```sh
git diff --stat 41054aa~1 e090909 -- testing.md parallel_refinement.md \
    project_structure.md docs/increments/03-pslg.md docs/increments/01-predicates.md
```

prints two files — `docs/increments/01-predicates.md` and
`parallel_refinement.md`. `testing.md` (items 8, 9, 10),
`docs/increments/03-pslg.md` (items 3, 4) and `project_structure.md` (items 5,
6, 7) were listed as fixed by that PR and were not touched by it.

**And the file-level diff overstates even the two that were touched**, which is
this section's own instance of the defect it is reporting. An earlier revision of
this paragraph read the two changed filenames as "items 1, 2, 2a, 2b, 2c and 11
shipped". Checked line by line rather than file by file, `parallel_refinement.md`
received **2a and 2b only**: step 6 still emitted "a list of `(p0, p1, is_river)`
segments" (item 1) and step 4 still said "at all intersections lying on it"
(item 2c), the exact exact-incidence phrasing the hot-pixel finding exists to
kill. A file appearing in a diff is not an item being fixed, in precisely the way
a cited grep is not a run one. Both are fixed here, with item 2 folded into the
rewrite of "Edge metadata".

All of it is re-found during this increment and is therefore this PR's, not a
ledger entry. `05-noder.md`'s own section now says so, at its head, and its stale
line numbers are marked stale rather than silently repaired.

**`testing.md`** — item 8, the `noding` invariant list, all five bullets, in the
replacement text `05-noder.md` specifies; item 10, line 211's
`noding_generators.h` marked as arriving at 5b, and the layout block's file list
naming `prop_noding_snap_invariants.cpp` (5a, shipped) beside
`prop_noding_no_crossings.cpp` (5b). Item 9 — the section's `[planned]` marker
becoming `[live]` — is **5c's**, not this PR's: the marker describes what CI
enforces over a module that ships, and 5b ships a noder that nothing calls.

**`docs/increments/03-pslg.md`** — item 3, lines 67-69's "sparse override set";
item 4, the `NodedPslg` promise at 486-493, which lists two guarantees and needs
the three it is missing, including that guarantee 9 does not survive. And one
`05-noder.md` did not catch: the same passage calls the wrapper signature change
"mechanical", which risk 13 refutes.

**`project_structure.md`** — items 5, 6, 7: the `noding` section's repeat of the
sparse-override-set claim; the directory listing, which has no
`include/terrain/noding/` and no `core/snap_grid.hpp` or `core/noded_pslg.hpp`,
and marks `src/noding/` *(planned)* when it is not needed at all; and the
broad-phase independence number, 4e12 buckets over a 100 km domain at 5 cm.

**`parallel_refinement.md`** — the "Edge metadata" section, whose claim that
non-river constraint edges need not remember their feature is false for every
linear feature that is not a river. Rewritten over property sets, with the
face-based half kept because it is sound for area features, and with
"geometrically forgotten" kept because it is true. This is the user's
correction, arriving mid-design; it is fixed where it lives rather than
recorded.

**`docs/increments/04-cdt.md`** — the paragraph beneath the mapping table saying
what stays reachable from a `NodedPslg`. **5c's**, for the same reason as
`testing.md`'s marker: it is false until the signature changes.

**`docs/increments/05-noder.md`** — five corrections, the first three of which
this document's premises establish, the fourth of which is the user's, and the
fifth of which came out of `@tester`'s red round:

1. Guarantee 14(a) is false on duplicate edges; the amended clause replaces it.
2. The 5b row of "Files and LOC" names no Python file and estimates the signature
   change at ~10; `bindings/core.cpp:466` refutes it.
3. The "Contingency split of 5b" section proposes a seam that this document does
   not take; it is replaced by a pointer here rather than left as a live
   alternative a later reader might act on.
4. The `is_river` representation paragraph is superseded in its *width* — a
   property set merged by union, not one bit merged by OR. Its reasoning about
   density, about the merge being keyed on node ids and about there being no
   single source bit to override survives verbatim, so it is marked rather than
   rewritten.

   *Historical, and true when written: the marking ends by saying this file's
   "Edge properties" section "carries the ruling and the scope". It carried both
   until increment 7 shipped the type; as of `6f7f7c8`, where this design was
   rebased onto that increment, it carries the scope and the shape on
   `NodedPslg`, and `docs/increments/07-edge-properties.md` is the ruling on the
   type. The clause is left standing as the record of what was true then;
   `grep -n 'carries the ruling and the scope' docs/increments/05-noder.md`
   finds it, and it belongs to whoever next edits that file.*
5. **Risk 3 names the wrong mitigation.** It says the winding re-check is what
   stands between a reversed sliver hole and a silent wrong mesh; measured, the
   refusal comes one stage earlier from `NonSimpleRing`. Corrected in place, with
   the mechanism and the residual pointing here.

**`docs/increments/03-pslg.md`** and **`project_structure.md`** each carry the
same widening where they describe the per-edge array, for the same reason.
