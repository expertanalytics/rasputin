# Increment 118: every --seed names its coordinate system

**Status: built** (as-built record, written by `@architect` on 2026-10-31
after the code; branch `worktree-seed-crs`, PR #413: `779eb03c` the change
and its tests, `fbaaf05c` the Bygdin script). Full lane by path
(`src_python/tin_engine/catchment.py`), Python only: no C++, no `@perf` run.
Number 118: 116 is the highest on master, 117 is taken on another branch
(`git ls-tree` over every ref and worktree).

## 0. The plain picture

- **Before:** `rasputin catchment` and `rasputin network-catchment` read
  `--seed X Y` as longitude and latitude (EPSG:4326) unless `--seed-crs`
  said otherwise. A seed given in metres without `--seed-crs` was read as
  degrees.
- **After:** `--seed` needs `--seed-crs`. A run without it is refused, and
  nothing is guessed. This is the rule increment 114 applied to `--bbox`:
  no coordinate's CRS is left implied.

      rasputin catchment --dem dtm10 --seed 8.5425 61.3512 --seed-crs EPSG:4326 --out c.geojson

## 1. Prior art: legacy and literature

*Legacy.* `git grep -n -i seed legacy-archive -- legacy` returns nothing.
Nothing to carry.

*Literature.* As for `--bbox` (increment 114, section 1): OGC API -
Features and WFS read a point or box with no CRS as WGS 84 longitude and
latitude. We depart from that default for the same reason: a missing CRS
is refused, not assumed. Axis order stays x then y (`always_xy`). No
novelty claimed.

## 2. What changed

- `CatchmentRequest.seed_crs` (`catchment.py`) has no default; building a
  request without it is a validation error. Its other callers already pass
  it (`catchment_batch.py`, the two in `cli.py`).
- `cli.py`: one shared `SEED_CRS` option (default `None`) for both
  commands, and one refusal, `NO_SEED_CRS`, worded as `NO_BOX_CRS`:
  "--seed needs --seed-crs, the CRS its two numbers are in (EPSG:4326 for
  longitude and latitude in degrees); a seed's CRS is never assumed".
- `network-catchment --unit X --seed-crs C` is refused with "needs --seed",
  like the other seed-only options there.
- The MCP tools follow without code: their arguments are derived from the
  click options (`mcp/schema.py`), so `seed_crs` becomes optional there and
  the command refuses it missing.
- `docs/benchmarks/2026-09-29/bygdin/run.sh` names `--seed-crs EPSG:4326`.
  The other benchmark script with a seed already named its CRS.

## 3. How it was tested

- `test_catchment.py`: a request without `seed_crs` raises a
  `ValidationError` naming the field; the existing requests name
  EPSG:4326.
- `test_cli_catchment.py`: a seed without its CRS is refused with the exact
  `NO_SEED_CRS` words; the seeds in degrees name EPSG:4326 and still give
  the same polygon, bit for bit, as the seed in the DEM's CRS.
- `test_cli_network_catchment.py`: the point seed names EPSG:4326;
  `--seed-crs` with `--unit` is refused.

## 4. Size

13 counted lines (`python3 tools/count_loc.py 571c6924 fbaaf05c`: 21 added,
8 removed). Old records that show `--seed` without `--seed-crs`
(increments 15c, 85, 106) are history and stay as written.
