# Increment 116: feature layers for plugins

**Status: designed** (`@architect`, 2026-10-30, on master `571c6924`,
branch `worktree-feature-table`). Full lane, one PR, Python only: no C++, no
`@perf` run. Number 116: no file under `docs/increments/` with that number
on any worktree or ref (`git worktree list`, `git for-each-ref` with `git
ls-tree`).

## 0. The plain picture

- **Today** a plugin can read a DEM through `tin_engine.plugin_api`, but not
  a feature layer with its attribute columns, and cannot write one.
- **After this increment** a plugin reads a whole GeoPackage or GeoJSON
  feature layer, every attribute column and its type, optionally only the
  rows near a box, and writes one with the core's own GeoPackage writer. It
  also finds a catalogue source's file the way `mesh` does.

## 1. Prior art: legacy and literature

*Literature.* No published method; two file standards. The GeoPackage
column types are those of OGC 12-128r18 (GeoPackage Encoding Standard
1.3.1), Table 1 under Requirement 5: `BOOLEAN`, `TINYINT`, `SMALLINT`,
`MEDIUMINT`, `INT`/`INTEGER`, `FLOAT`, `DOUBLE`/`REAL`, `TEXT{(n)}`,
`BLOB{(n)}`, `DATE`, `DATETIME`. All but `BLOB` are mapped below; `BLOB` is
refused. GeoJSON is RFC 7946: `properties` is any JSON object and a
feature's `id` member is optional (section 3.2), so field types are
inferred from the values and the fid is the feature's position. Field type
names (`esriFieldTypeInteger` and kin) are the ones the core's writer
already takes (`src_python/tin_engine/io/geopackage_write.py@571c6924:40`).
No novelty claimed; none searched.

*Legacy.* `git grep -l -i -E "geojson|gpkg|geopackage|shapefile|ogr|fiona"
legacy-archive -- legacy` lists five files; four match only "geographic"
through `ogr`. The fifth, `legacy/rasputin/gml_repository.py`, parses one fixed GML schema with
`lxml` into a hard-coded class list, without a generic column read. Nothing
to carry.

*Ideas file.* No entry taken up or overlapped.

## 2. Design

A plugin may read a whole feature layer, with every attribute column, and
write one, through `tin_engine.plugin_api`.

```python
# tin_engine/io/feature_table.py (new)
@dataclass(frozen=True, slots=True)
class FeatureTable:
    """A feature layer as read: its CRS text, its GeoPackage layer (None for
    GeoJSON), its columns as (name, field type) in file order, and per row
    (fid, geometry or None, values in column order), in file order."""
    crs: str
    layer: str | None
    columns: tuple[tuple[str, str], ...]
    rows: tuple[tuple[int, BaseGeometry | None, tuple[Any, ...]], ...]

def read_feature_table(
    path: Path,
    layer: str | None,
    box_for: Callable[[str], tuple[float, float, float, float]] | None = None,
) -> FeatureTable: ...
```

- A GeoPackage: the layer (`layer_info`), its columns less key and geometry
  (`layer_summary`), every row or, with `box_for(own CRS)`, the rows its
  R-tree query returns (`query_features`, all columns as `extra`). Declared
  type to field type: `INTEGER`, `INT`, `MEDIUMINT`, `SMALLINT`, `TINYINT`,
  `BOOLEAN` to `esriFieldTypeInteger`; `REAL`, `DOUBLE`, `FLOAT` to
  `esriFieldTypeDouble`; `TEXT`, `TEXT(n)`, `DATE`, `DATETIME` to
  `esriFieldTypeString`; any other is a `FeatureError` naming the column.
- GeoJSON (`.geojson`, `.json`): the properties' keys in order of first
  appearance; a key's non-null values all `int` (not `bool`) to Integer, all
  `int` or `float` to Double, all `str` to String; `bool` values to Integer;
  any other mix is a `FeatureError` naming the key. The fid is the 1-based
  position; a feature's own `id` member is not kept. `box_for` is not used
  (no index).
- Any other suffix is a `FeatureError`, as `read_source`'s.
- The writer is the core's existing one, re-exported: `srs_entry`,
  `create_geopackage`, `create_features_table`, `insert_features`, which
  `gml-to-gpkg` uses today.
- `sources.features_file(source, cache) -> Path`: a features source's own
  `path`, else `<cache>/<id>/features.gpkg`; a `ValueError` naming the
  missing copy, or the missing cache. `cli._cached_features` calls it and
  keeps its own refusal words (with the fetch to run). `cache_root` moves
  from `cli` to `sources` unchanged, and `cli` imports it from there.
- `plugin_api` re-exports `FeatureTable`, `read_feature_table`,
  `FeatureError`, `srs_entry`, `create_geopackage`, `create_features_table`,
  `insert_features`, `catalogue`, `user_source_file`, `plugin_source_files`,
  `RemoteSource`, `cache_root`, `features_file`, `kept_named` (the rule of
  an entry's `named_by`, `sources.kept_named`, unchanged); `docs/plugins.md`'s
  table gains a row for each.

Blocking (sqlite3, GEOS): an async caller runs the reader in
`asyncio.to_thread`. No behaviour of any core command changes.

No new constant.

## 3. Tests (`tests/python/`)

Each fixture is built in the test; its figures are the ones it is written
with.

- **P1.** A GeoPackage written with `create_features_table` /
  `insert_features` in EPSG:25833, fids 7, 9, 12, columns `code` String, `n`
  Integer, `x` Double: `read_feature_table` returns those columns and types,
  those fids, the values exactly, each geometry `equals_exact` at 0, and the
  CRS text.
- **P2.** The same file with a `box_for` meeting the boxes of rows 7 and 9
  only: two rows.
- **P3.** GeoJSON `{"kind": "a"}` then `{"kind": "b", "n": 2}`: columns
  `kind` String, `n` Integer; row 1's `n` is None; fids 1, 2. A key that is
  `"1"` in one feature and `1` in another: a `FeatureError` naming it.
- **P4.** A GeoPackage with a BLOB column: a `FeatureError` naming it.
- **P5.** `features_file`: a local entry's path; a remote features entry
  under the cache; a missing copy and a missing cache, each a `ValueError`;
  the existing `mesh` refusal tests pass unchanged.
- **P6.** Each new `plugin_api` name is importable from it and named in
  `docs/plugins.md`.

No suite here is invariant-critical; no mutation round.

## 4. Size

About 70 counted lines (estimate): `feature_table.py` about 50,
`features_file` and the `cache_root` move about 5 net, `plugin_api` about
16.
