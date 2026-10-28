# Increment 110: the preferred source's shorelines are the land cover's lines

**Status: designed** by `@architect`, 2026-10-28, on master `aabac548`,
worktree `direct-lake-density`. One PR, full loop (about 85 net production
lines), Python only. Number 110: 109 is the highest under `docs/increments/`
over every ref.

Ola, 2026-10-28: "We're still not quite there on the direct MCP, two stage
looks much better :)"

## 0. The plain picture

The cause is in `docs/research/sula-direct-lake-shore-density.md`: with
`--prefer lake=nve-lakes`, a direct run takes the lake *classes* from NVE but
keeps CORINE's lake *outlines* among the constraint lines, beside NVE's
shorelines. Each lake gets two crossing shorelines, and the start-quality step
fills the strips between them with small triangles. After this increment,
when the `--land-cover` file is also a `--features` source, its polygons get
the preferred classes **before** they become lines, by the same
`preferred.overlay` rule the class fractions use (107, section 2). Each lake
then has one shoreline, NVE's, simplified together with its neighbours in one
coverage, and the lines bound exactly the classes the fractions assign. The
direct command becomes:

    rasputin mesh ... --features <clc2018 features.gpkg> --features-rules <corine lines>
        --land-cover clc2018 --prefer lake=nve-lakes

No new option. A run without `--prefer`, or one whose land cover is not a
`--features` source, is unchanged.

## 1. Prior art: legacy and literature

*Legacy.* `git grep -n -i -E "shore|constraint.*lake|lake.*constraint" legacy-archive -- legacy`:
no files. One source per mesh there; nothing to carry across.

*Literature.* Replacing one map's features of a theme with a better source's
and repairing the topology around them is map conflation (Saalfeld,
"Conflation: automated map compilation", *Int. J. Geographical Information
Systems* 2(3), 1988, 217–228). Here the replacement rule is 107's, and the
shared borders are then simplified once, as one polygon coverage (the GEOS
coverage simplifier the clean-up already uses), so two sources cannot give
two near-copies of one border. No novelty is claimed; none searched.

*Ideas file.* Takes up no entry of `../rasputin_scratch/ideas/ideas.md`.

## 2. Blueprint

```
CLI: --land-cover file == a --features file (resolved path, same layer)?
   │ yes, and --prefer given
   ├─ read_polygons(each --prefer source, DEM CRS, domain outline)   (as 107 reads them)
   └─ open_features(request, domain, dem_crs, cell=..., prefer=Prefer(...))
         _Tally._take, that source only: polygons collected as cover, moved to DEM CRS
           → preferred.replace(base, chosen, near=[read region])  ← the cut and refill of 107
           → _clean (repair, same-class merge, simplification, as asked)  → lines
```

- `tin_engine/preferred.py`: split `overlay` into
  `replace(base, chosen, near: Sequence[BaseGeometry]) -> Replaced` (steps
  1–4 of 107: keep what meets `near`, cut, leftover pieces, neighbour; a
  frozen dataclass of the kept, given and refilled polygons, with
  `polygons`) and `overlay`, which calls `replace` with the triangles and adds
  the areas. `overlay`'s results do not change. About 15 lines.
- `tin_engine/feature_input.py`: a frozen dataclass
  `Prefer(source: int, rules: MeshingRules, chosen: Sequence[tuple[int, Sequence[BaseGeometry]]])`:
  the index of the land cover's source in `request.sources`, the land cover's
  meshing rules (value to class code), and 107's `chosen` in the DEM's CRS.
  Data, so it is an argument of `open_features`, not a field of the Pydantic
  `FeatureRequest`. In `_take`, for that source: every polygon goes the
  clean-up way even when every clean-up switch is off (the stage then only
  clips and runs `coverage_clean` at 0); each polygon's overlay code is its
  value read by `_code` with the land cover's rules; after the move, before
  `_clean`, the cover is `replace`d with `near` = the read region. Each output
  polygon takes the mask, code and merge key of its class in that source, a
  given polygon those of the replaced class (`_mask` of the code's text), and
  the fid of its class's first polygon (None where the source has none). About
  40 lines.
- `cli.py`: find the land cover's `--features` source; read the preferred
  polygons with the domain outline's points; pass `Prefer`; the refusals
  below; one sentence in `--prefer`'s help ("When the --land-cover file is
  also a --features source, its lines take the same classes."). About 30
  lines.

**Refusals** (usage errors, before any phase is timed):
1. A `--features` file that is also a `--prefer` source, while the land
   cover is a `--features` source: its shorelines already come in through the
   land cover's lines, and a second copy is the defect this fixes. Names
   `--features`, the file and the key.
2. The land cover's `--features` rules read another attribute than its
   `--land-cover` rules (the codes would not match), or keep the polygons out
   of the clean-up (`cleaned` false). Names `--features-rules`.

## 3. Tests

1. `replace` with `near` = one square gives `overlay`'s polygons for 107's
   test 1 (base 100×100 m class 322, 512 square [40,60]², lake
   [30,50]×[40,60]); `test_preferred.py` stays green unchanged.
2. `open_features` on that square, domain = the square, clean-up off and on
   (repair 0, merge on, tolerance 0): the lines off the outline are the lake
   rectangle's boundary, 80 m, every vertex within 1e-9 m of it; none of the
   60 m of the old 512 border off that boundary (x = 40 inside the lake, and
   the east part beyond x = 50) remains.
3. CLI on 107's scene (`test_cli_mesh_prefer.py`: `cover.geojson` as
   `--land-cover` with `corine` and as `--features` with a non-coded
   lines file, `lake.geojson` preferred, `--features-tolerance 0
   --features-repair 0 --features-outline-snap 0`). Measured with
   `preferred.overlay` and shapely: the lines inside the domain total
   9,838.4 m today and 9,082.7 m after; the lake's 2,600.4 m border inside the
   domain is covered (0.1 m of it today); none of the 512 polygon's 2,756.4 m
   border inside the domain is kept; 3,350.1 m of today's lines lie farther
   than 1 m from the new ones. The cell shares of 107's case 7 are unchanged.
4. Refusal 1 (`--features lake.geojson` added) and refusal 2 (other
   attribute), exit 2, naming the option.
5. `--prefer` with only `lake.geojson` as `--features` (land cover not among
   them): the lines layer equals the run without `--prefer`.

## 4. Acceptance on Sula (NVE 73.27.0), 20 m

The direct command of `docs/research/sula-direct-lake-shore-density.md` less
`--features nve-lakes` and its rules. Emulated without code by
`docs/increments/110-probes/direct_prefer_probe.py` (CORINE overlaid with NVE
lakes written to a GeoJSON, meshed as the only `--features` source; the probe
reproduces today's run exactly):

| Run | Triangles | Small, in band | Small, elsewhere | Area km² | Lake |
|---|---|---|---|---|---|
| direct, today | 6,429 | 1,898 | 125 | 30.4411 | 8.935 % |
| direct, fix emulated | 3,271 | 184 | 94 | 30.4411 | 8.935 % |
| two-stage, same 15° start angle | 3,564 | 213 | 248 | | |

Small = under 300 m²; band = within 40 m of both a CORINE border and an NVE
shoreline (1.254 km²). **Accepted when** triangles ≤ 3,564 and small in band
≤ 213 (the two-stage run's, at the same start angle; about 3,271 and 184
expected: the real route clips and repairs slightly differently from the
emulation), and area 30.4411 km² and lake 8.935 % as printed to those digits.
Then the same through the MCP server's mesh tool.

## 5. Size, timing, what is left out

About 85 net production lines. Added time, only on `--prefer` runs whose land
cover is a `--features` source: a second read of the preferred source and a
second cut and neighbour search (Sula: 0.03 s and 0.13 s, measured; on
107's Numedalslågen figures about 2.9 + 6.4 s). No `@perf` acceptance run:
no refine or mesh code changes, and `tools/bench.py`'s cases pass no
`--prefer`, so their runs are the same by construction. Left out: reading the
preferred source once for both uses, and `--prefer` changing a `--features`
source other than the land cover's own file.
