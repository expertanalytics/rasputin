# Increment 9 — gallery output

Status: design. No production code and no test on this commit.

Depends on increment 8 (`docs/increments/08-crossing-gallery.md`), which is
merged into this branch and adds `road-enters-forest`, `wall-leaves-domain` and
`bridge-over-lake`. Depends on no C++ change and adds no binding.

## The scope was cut, and this says why

The first version of this record was four times the size. It committed eleven
rendered SVGs under `docs/gallery/`, added an HTML index page, added a
byte-comparison test to stop the committed files drifting from the code, and
argued at length with `06-cdt-viewer.md` about whether that comparison was the
golden-file test that record rejects.

The user cut all of it:

> *"Don't waste a lot of energy and tokens on svg-visualization. There are
> plenty of good tools to render TINs, and standard formats to save triangular
> meshes… Simple vis of smaller concepts using svg is fine, at the moment. But
> for surfaces, this approach is silly."*

That is the rule to keep. **SVG here is for small concepts, not for surfaces.**
The gallery is eleven hand-typed shapes of a dozen vertices each, and drawing
those as SVG is right. Drawing a real terrain surface that way is not, and
neither is building a publishing pipeline — an index page, committed artifacts,
a drift check — around pictures of test shapes.

So do not reinstate any of it. If a surface needs looking at, the answer is a
mesh file and a tool that reads meshes; `docs/increments/10-mesh-output.md` is
that answer. The gallery's job is smaller: render the eleven fixtures where the
caller asks, in one command.

`06-cdt-viewer.md`'s outright rejection of a golden-file SVG comparison stands
as written. Nothing here amends it.

## The defect

In the user's words: *"the rasputin-gallery, which quite surprisingly landed
files outside the repo, or temporary location, which is not very impressive,
has not been updated with new examples"*.

Two facts, each checked rather than asserted.

**Every picture lands outside the repository.** `cli._destination` with no
`--out` returns a path under a fresh `tempfile.mkdtemp(prefix="rasputin-")`.
The comment above that line states the reason, and the reason is sound for what
it covers. It is the *only* mode that exists, and that is half the defect.

**There is no way to render the set.** `draw` takes one fixture name. Eleven
fixtures is eleven commands and eleven temp directories.

```sh
.venv/bin/python -c "from tin_engine.viz.fixtures import GALLERY; print(len(GALLERY))"
```

## Prior art in `legacy/`

```sh
grep -rlniE "write_scene|\.html|data\.js|svg|gallery|thumbnail" legacy/
```

returns two files:

```
legacy/rasputin/geometry.py
legacy/rasputin/web_visualize.py
```

`geometry.write_scene` copies a bundled web-template directory to an output
directory and writes `data.js` beside it. `web_visualize.visualize_tin` calls
it and prints instructions to start `python -m http.server 8080`.

**Nothing is carried across.** One thing is carried across as a negative:
`write_scene` begins `if output.exists(): shutil.rmtree(output)`. A command
that deletes a directory the user named is the opposite of this project's
output boundary. No command designed here gains the power to delete.

The legacy renderer is WebGL over a real TIN. It answers the surface question,
not this one, and it is the reason the surface question is increment 10's.

## Ruling 1 — a second command, not a flag on `draw`

**`rasputin gallery` renders every fixture into one directory. `draw` is
unchanged.**

`--all` was the cheaper-looking option and it is wrong, because it makes
`--out` mean two things: a file with a fixture name, a directory without one.
An option whose type depends on another option has to be validated at runtime
with a message explaining which mode the user is in.

Surface:

```
rasputin gallery [--out-dir DIR] [--out-parent DIR]
                 [--delaunay/--no-delaunay] [--snap-spacing X]
```

`draw`'s per-picture knobs — `--labels`, `--label-limit`, `--vertices`,
`--title` — are not repeated here. They inspect one picture; this command
produces a set. `--delaunay` and `--snap-spacing` stay, because they are the
two levers worth comparing across the whole set at once.

Filenames are `<fixture>.svg`, one per `GALLERY` key.

Exit code 0 whenever pictures were produced, exactly as `draw` does, including
for the two deliberate failure fixtures. A drawn failure is a drawn picture.
The command prints the directory last.

## Ruling 2 — the temp-directory default survives, on both commands

**`draw` with no `--out` keeps writing to a fresh temp directory. `gallery`
with no `--out-dir` does the same.**

The stated reason still holds: a second run must not overwrite a picture still
open in a browser. It was never wrong; what was wrong is that it was the only
mode. Ruling 1 supplies the other mode.

Making `--out-dir` mandatory was the alternative, and is rejected: the
exploratory run — "what does the whole set look like at this spacing" — is the
common one, and forcing it to name a directory makes a `--snap-spacing 1.0`
experiment land somewhere permanent.

The invariant, stated once:

> **No rasputin command writes inside the repository unless a path on the
> command line says so.**

## Ruling 3 — the directory boundary extends `_destination`, and weakens nothing

`cli._destination` is a hostile boundary and its docstring says which attacks
it turns away. `gallery` needs a directory rather than a file, so it gets a
sibling — `_destination_dir(out_dir, out_parent) -> Path` — applying the same
rules in the same order:

- `None` means a fresh `mkdtemp`, per ruling 2.
- A symlinked final component is **refused, not followed**. Same reason as
  `_destination`'s: a link pointing back inside a permitted parent passes
  containment and still writes where the user did not ask.
- `resolve()` first, then containment against `--out-parent`, so `..` cannot
  smuggle a write out of a permitted parent.
- A missing directory is a message, not a traceback, and the command does not
  create it.

Two rules are **added**, because a directory destination has two exposures a
single named file does not:

1. **Each filename is checked before it is joined.** A fixture name must match
   `^[a-z][a-z0-9-]*$`. Every `GALLERY` key satisfies this today, which is
   exactly why the check is free now — and the day someone adds a fixture
   called `../x` it is the difference between a message and a write.
2. **An existing entry that is a symlink is refused rather than overwritten.**
   `_destination` refuses a symlinked `--out`; without this, that protection is
   bypassed by writing into a directory instead of at a file.

Nothing is removed and nothing is loosened.

## The blueprint

Nothing new crosses any boundary. The command is a loop over what `draw`
already does:

```
GALLERY       -> cli._triangulated   -> Attempt      (the only step touching _core)
Attempt       -> viz.build_scene     -> Scene
Scene + style -> viz.render_svg      -> str          (11 of these)
11 strings    -> cli writes files                    (the only I/O)
```

`viz/` is untouched: no new module, no new import, and `06-cdt-viewer.md`'s
"no file is written below `cli.py`" stays literally true. `cli.py` stays the
single composition root.

## Files and LOC

**Estimate, not a measurement** (`PRINCIPLES.md` B4).

| File | What | Est. non-comment lines |
|---|---|---|
| `src_python/tin_engine/cli.py` | `gallery` command | ~18 |
| `src_python/tin_engine/cli.py` | `_destination_dir` and the two added checks | ~14 |
| | **total** | **~32** |

About 5% of `CLAUDE.md` §2's ceiling. No seam to pre-declare, no new module.

Measure it with the instrument `06-cdt-viewer.md` and `08-crossing-gallery.md`
used, noting as those files do that it counts docstring lines and has a known
blind spot on a bare `*,`:

```sh
git diff master...HEAD -- src_python/tin_engine \
  | grep '^+' | grep -v '^+++' | sed 's/^+//' \
  | grep -vcE '^\s*(//|#|\*|/\*|\*/|$)'
```

Increment 8 is merged into this branch, so master is the right base.

Tests are excluded from the ceiling. Estimated at ~50 lines: `gallery` writes
one file per fixture and prints the directory, and `_destination_dir`'s four
refusals.

## Testing

No invariant-critical suite, so no mutation round, per the cost constraint in
`docs/increments/README.md`. There is no topology decision here.

`_destination_dir` gets the same four negatives `_destination` has: a symlinked
directory, a directory outside `--out-parent`, a missing directory, and — new
here — a symlinked entry inside an otherwise fine directory. A fifth covers the
filename pattern, with a name the pattern rejects.

The `gallery` command's test writes into `tmp_path` and asserts one file per
`GALLERY` key, including the two failure fixtures.

## Acceptance

A person runs

```sh
rasputin gallery --out-dir /some/dir
```

and gets eleven SVG files in that directory, one per fixture, failures included.

## Not in scope

- **Committed pictures under `docs/gallery/`, an index page, and any check that
  a committed picture still matches the code.** Cut by the user; the reasoning
  is at the top of this record. Do not reinstate.
- **PNG or any raster output.**
- **A served page, a viewer, or any JavaScript.**
- **Anything about surfaces.** `docs/increments/10-mesh-output.md` holds that.
- **Any change to `_core`, the bindings, the stub file, or `viz/`.**
- **Rendering anything but the gallery fixtures.** Real data has no ingestion
  path yet; `ROADMAP.md`'s MVP-gap list holds that.

## The record this updates

`ROADMAP.md` gains a row for increment 9, in this PR, per
`docs/increments/README.md`'s merge rule.

`06-cdt-viewer.md` gains a short pointer to this file, **appended at the end of
that document rather than inserted into its status block.** That is deliberate:
`07-edge-properties.md` cites `06-cdt-viewer.md` by line number in three
places, and inserting near the top would leave all three resolving to the wrong
text while still resolving. Appending shifts nothing. Find the citations with:

```sh
grep -rn "06-cdt-viewer.md:" docs .claude src_python tools
```

Nothing in `06-cdt-viewer.md` is made false by this increment, and no ruling in
it is amended.
