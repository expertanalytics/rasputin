# Increment 78: MCP over HTTP, inputs by link and by upload

**Status: designed** by `@architect`, 2026-10-19, on `af9f5b52` (master,
increments 69 and 72 in), branch `worktree-mcp-inputs`. Full path (about 130
counted lines, over the short path's 50); no C++, no refine or mesh code, so
no `@perf` run. Number 78 because 70, 71 and 73 to 77 are held by designs on
other worktrees (`ls .claude/worktrees/*/docs/increments/`).

## 1. What Ola asked

As the main session reports Ola's yes of 2026-10-19 (no verbatim quote was
passed on): over HTTP (`rasputin mcp --http`, increment 69), a tool input may
be (a) one of the server's own download links, and (b) a file the client
uploads, so a client with no shell and no shared folder (Codex limited to
URIs) can run a whole exchange by URI: upload a polygon, get a link, mesh
with that link as `--domain`, get the mesh's link back. Today a link given
as an input is refused ("not a plain file path", `src_python/tin_engine/mcp/workspace.py@af9f5b52:73`),
inputs are paths in the user's `work/`, and nothing can be uploaded.

## 2. Prior art: legacy and literature

*Literature.* MCP revision 2026-07-28 has no file input: tool arguments are
JSON. The File Uploads Working Group (charter of 2026-04-23,
modelcontextprotocol.io/community/file-uploads/charter, read 2026-10-19)
says "servers that need a file from the user resort to prose instructions
asking for base64 strings or local paths"; its first proposal, SEP-2356
(RFC 2397 data URIs as tool arguments), was closed unmerged on 2026-06-26 in
favour of a successor not yet accepted. **Departure:** with no standard, this
step takes the content as a plain tool argument (text, or base64 for
binary), the practice the charter describes; when the working group's form
lands, `upload` can accept it too. A link given back as an input follows the
spec's non-normative *Stateful Tools* section, already cited by increment
69: "For authenticated servers, a handle is a name, not a capability. The
server should validate the caller's authorization against the handle on every
call." So a link as an input is judged by the caller's bearer token and the
user in its path, not by its signature or expiry (section 5.1). Novelty:
none claimed.

*Legacy.* Nothing carries across:

    $ git grep -l -i -E "upload|b64decode|base64|data:[a-z]|urlparse|urlsplit" legacy-archive -- legacy
    (no output, exit 1)

*In this tree.* Increment 69's `CallScope`, `refusal`, `Store`, `Download`,
`_download_link` and the `NAME` and `USER` patterns of `http.py`; increment
72's `argv` (positionals after `--`) and the sidecars set through
`model_copy(update=...)`, the pattern section 5.1 reuses.

## 3. Measured

- **Input sizes** (`find ../rasputin_scratch ../rasputin_data -name '*.geojson'`,
  `stat -f %z`, 2026-10-19): 5,986 files, median 36,282 bytes, 99th
  percentile 934,228; 5,945 are at most 2 MiB. The 41 larger ones are
  land-cover extracts and summaries of 55 to 74 MB, data that belongs among
  the server's sources, not in an upload.
- **The request limit:** the SDK answers a body over 4 MiB with 413 before
  any tool runs (`DEFAULT_MAX_REQUEST_BODY_SIZE`, `mcp/server/transport_security.py`
  line 15 and lines 147 and 163, `mcp` 2.3.0). Base64 of 2 MiB is 2,796,204
  characters, so a 2 MiB upload fits with room for the JSON-RPC envelope.
- **Input formats name no other file:** GeoTIFF is read by tifffile alone
  (side files ignored, `src_python/tin_engine/io/repository.py@af9f5b52:138`),
  a GeoPackage is opened `mode=ro` (`src_python/tin_engine/io/repository.py@af9f5b52:514`), legends
  and palettes are TOML with no file keys (`grep -n "Path(" src_python/tin_engine/legends.py src_python/tin_engine/palettes.py`
  finds only the value given). So an uploaded file cannot make the child open
  another (increment 69's round-3 review asked for this check).

## 4. Blueprint

```
call_tool(ctx, params)                              [server.py]
  validate → scope = scope_for(ctx)
  _inputs(scope, tool, arguments)                   NEW, before both checks
     each write value that is a URI      → refusal (HTTP only)
     each read value: workspace.linked(scope, argument, value)
        None (not a URI, or stdio)       → unchanged
        refusal words                    → refused, child not run
        Linked(text="../store/<name>", place)
                                         → place must be a regular file (lstat)
                                         → os.utime(place): kept 7 more days
     tool.substitute(arguments, {value: text})       [schema.py, pure]
  check() → slot → check() → child, as today         (unchanged)
upload(name, text | base64)                         NEW tool, HTTP only
  name, size, encoding → admit (sweep, quota)
  → Store.receive(name, data): work/uploads/<name>  [store.py]
  → delivery.deliver → one Download + resource_link (as download_link)
```

Module rules, beyond increment 69's:

- **`workspace.py`**: `CallScope` gains `link_prefix: str | None = None`,
  `public_url.rstrip('/') + '/files/'`, as `links.url` builds it, over HTTP, None on stdio. New pure
  `linked(scope, argument, value) -> Linked | str | None` and the frozen
  dataclass `Linked(text: str, place: Path)`. No I/O.
- **`schema.py`**: `DerivedCommand.substitute(arguments, mapping) ->
  BaseModel`: the same walk as `reads` (each item of a repeatable value,
  each string of a tuple), set fields only, every value equal to a key
  replaced; validated again from `model_dump(exclude_unset=True)`, so
  `model_fields_set` and the argv order are unchanged.
- **`store.py`**: `Store.receive(name, data: bytes) -> Path`: writes
  `work/uploads/.part-<hex>` then `os.replace` onto `work/uploads/<name>`,
  in a worker thread; refuses (`ValueError`) when `work/uploads` resolves
  anywhere but itself.
- **`server.py`**: `_inputs`, `UPLOAD_TOOL`, `_upload` (beside
  `_download_link`, sharing its result shaping), the call log line for an
  upload.
- **`http.py`**: `scope_for` sets `link_prefix` and adds `<user>/store` to
  `read_roots`; `build_app` refuses to start when a command is named
  `upload`; `HTTP_INSTRUCTIONS` gains two sentences (section 5.3).

## 5. Decisions

### 5.1 A link as an input

A read value is **a URI** when it contains `://` or starts with `data:`
(any case). Over HTTP, every URI is judged by `linked`, with no I/O:

1. It must start with `link_prefix` exactly, as a string (no URL parsing,
   no normalisation, case-sensitive). Otherwise refused: *not a link this
   server made*. So `https://mcp.example.org.evil/files/...`,
   `https://mcp.example.org/filesx/...`, `file:///etc/hosts` and `data:`
   values never reach the file system, and **the server fetches nothing**.
2. The rest, cut at the first `?` or `#`, must be `<user>/<name>` with
   `USER` and `NAME` (`http.py`'s patterns, moved to `workspace.py`) each
   matching in full. `NAME` is `[0-9a-f]{64}(\.[a-z0-9]{1,8})?` and has no
   `/`, `.` run or `%`, so `../x`, `%2e%2e` and `<sha>.gpkg/../x` fail.
   Refused: *not a link this server made*.
3. `<user>` must be `scope.user`. Refused: *another user's file*.
4. The query (`fn`, `exp`, `sig`) is ignored: the bearer token has already
   said who calls, and the file is the caller's own (section 2). An expired
   link, and a link from before a restart, therefore work as inputs while
   the stored file is kept.
5. Result: `Linked(text="../store/<name>", place=<out-root>/<user>/store/<name>)`.
   The child's cwd is `<user>/work`, so the relative text names the stored
   file, is always a plain path, and resolves inside the new read root
   `<user>/store`. `_inputs` then requires `place` to be a regular file by
   `lstat` (a link or a swept file is refused: *no longer kept*) and touches
   its mtime, so the sweep at `admit` cannot remove it before the child runs.

The stored file keeps its suffix, which is how the CLI picks the format.
`store/` is written only by the server (increment 69, 5.1), so no call can
change the bytes a link names. Over HTTP, a URI in a **write** value (`out="https://…"`)
is refused over HTTP: a link names a stored file, which no call writes.

On **stdio**, `link_prefix` is None and `linked` returns None for every
value: a link given there goes to the CLI as a string, as today. Stdio has
no links of its own and shares the client's disk, so nothing changes there.

### 5.2 Upload: the `upload` tool

Listed over HTTP only, beside `download_link`. Arguments: `name` (required),
and exactly one of `text` (written as UTF-8) or `base64` (standard alphabet,
padded, `base64.b64decode(..., validate=True)`).

- **Name:** `[A-Za-z0-9_-][A-Za-z0-9._-]{0,127}`, in full: no `/`, no
  leading dot, so neither `.`, `..` nor a hidden file. The suffix picks the
  format the CLI reads it as (`.geojson`, `.gpkg`, `.tif`, `.toml`).
- **Size:** at most `UPLOAD_BYTES` = 2 MiB (2,097,152) after decoding,
  checked before decoding (text: its UTF-8 length; base64: at most 2,796,204 characters) and again after decoding. Scale: a catchment
  polygon; checked against the 5,986 GeoJSON files of section 3. A body over
  4 MiB never reaches the tool (the SDK's 413).
- **Quota:** `delivery.admit` first, as for a command, so a folder at or
  over the quota takes no upload.
- **Where it lands:** `work/uploads/<name>`, replacing a same-named earlier
  upload whole, then stored as `store/<sha256><suffix>` by `deliver`. The
  result is one `Download` (`path` = `uploads/<name>`) and one
  `resource_link`, as from `download_link`. So the client may pass either
  the path or the link; the earlier upload's link still names its own bytes.
- **Lifetime:** as every file of the user's folder: 7 days from the last
  write or use (5.1 step 5); the link itself 1 hour for downloading, any age
  as an input.
- **Annotations:** not read-only, not destructive (earlier bytes stay
  reachable by their link), idempotent, not open-world.
- **Log:** `upload user=<u> file=<first 12 hex> bytes=<n>`; never the name's
  content, never the link.

### 5.3 Errors, in plain words

| case | words (`<arg>` is the argument, value echoed only where shown) |
| - | - |
| URI not this server's | `<arg>: <value> is not a link this server made; the server opens no other address. Send the file with the upload tool and pass the link it returns.` |
| another user's link | `<arg>: this link names another user's file; a call reads only your own.` |
| stored file gone | `<arg>: the file this link names is no longer kept (files are kept 7 days after their last use); upload it again, or rerun the call that made it.` |
| URI as an output | `<arg>: a link names a stored file, which no call writes; give a path in your folder.` |
| upload name | `name: give a file name of letters, digits, '.', '_' and '-', not starting with '.', up to 128 characters.` |
| neither or both | `give the file as text or as base64, one of the two.` |
| too large | `the file is <n> bytes, over the upload limit of 2097152 (2 MiB); a larger input must be one of the server's data sources.` |
| bad base64 | `base64: not valid base64 (A-Z, a-z, 0-9, + and /, padded with =).` |

`HTTP_INSTRUCTIONS` gains: "An input may also be a download link this server
gave you, of any age, while its file is kept. upload sends a file of your
own, up to 2 MiB, and returns its link."

### 5.4 Codex's timeout (documentation only)

Codex stops waiting for a tool after 60 s by default (its
`tool_timeout_sec`), which cuts off a long mesh. The README's "Serving over
HTTP" part already says so and sets `tool_timeout_sec = 3600` (`README.md@af9f5b52:275`,
from increment 72), so this step only adds one paragraph there: links and
uploads as inputs, the 2 MiB limit.

## 6. Tests (`@tester`), red first

In `tests/python/test_mcp_http_inputs.py`, `pytest.importorskip("mcp")`,
through `mcp_fake_child.py` (its `argv_file` and `copy` keys) and
`mcp_http_driver.py`. No kernel except T7. **Invariant-critical: T1 and
T3**, mutation targets: `fullmatch` to `match` for `NAME` (killed by
`<sha>.gpkg/../x`, `<sha>0` (65 hex) and `<sha>.gpkg%2f..%2fx`), the user check dropped (Kari's case), the prefix test
replaced by `"/files/" in value` (the look-alike host), the `lstat`
regular-file test dropped (a symlink planted in `store/`), the leading-dot
rule dropped (`..`).

- **T1 `linked`, pure** (`public_url` `https://mcp.example.org`, user `ola`,
  sha of 64 hex). Accepted as `../store/<sha>.gpkg`: Ola's link; the same
  with `exp` in the past; with a changed `sig`; with no query. Refused as not
  this server's, each without touching the disk (the out-root does not
  exist): a host of the same length, `https://mcp.example.com/files/ola/<sha>.gpkg`, `https://mcp.example.org.evil/files/...`,
  `https://mcp.example.org/filesx/ola/<sha>.gpkg`, `http://` for `https://`,
  `file:///etc/hosts`, `data:text/plain,x`, names `../x`, `%2e%2e`,
  uppercase hex, `<sha>.gpkg/../x`. Refused as another user's: Kari's link.
  None: `glo30`, `EPSG:3006`, `in.geojson`. On a stdio scope: None for
  every value above.
- **T2 `substitute`, pure:** a link replaced in a plain option, an item of
  repeatable `features`, the path of `tolerance_near` (`(path, float)`) and
  a positional; `model_fields_set` and `argv` otherwise unchanged.
- **T3 a link round trip over HTTP.** A fake-child call writes `a.geojson`
  (1,024 bytes); its link, passed as `domain` to a second call, reaches the
  child as `../store/<sha>.geojson` (`argv_file`), and the child's `copy`
  of it to `b.geojson` holds the same bytes. The stored file's mtime is
  later than before. Then, each with the child not called: the stored file
  removed by hand (*no longer kept*); replaced by a symlink to a file in
  `work/` (*no longer kept*); Kari's call with Ola's link; a second app on
  the same out-root (new key) accepts the old link; `out` set to a link.
- **T4 upload.** Text, 1,024 bytes of GeoJSON: `Download.path` is
  `uploads/x.geojson`, its sha256 is the bytes', the file is in
  `work/uploads/` and `store/`, the link downloads the same bytes. Base64
  of bytes 0 to 255: same round trip. Exactly 2,097,152 bytes of text
  accepted; 2,097,153 refused, nothing written. Base64 of 2,097,153 bytes refused; of 2,097,152 accepted. Refused, nothing written:
  both or neither content, names `..`, `.x`, `a/b`, `a b`, 129 characters,
  base64 `@@@`, a folder at the quota (`--quota 1M`, 1 MiB in `work/`).
  The same name twice with other bytes: `work/uploads/` holds the second,
  and each link downloads its own.
- **T5 stdio unchanged.** No `upload` listed; a link given as `domain`
  reaches the child unchanged; `out` set to a URI is not refused by the new check; increment 37's and 69's MCP suites pass
  without edits.
- **T6 start and log.** A plugin command named `upload` stops start. With
  `caplog` on `tin_engine.mcp.http`: an upload logs one line with the user
  and no record holds the uploaded text or a `sig`.
- **T7 the whole exchange, real `mesh`.** Upload a polygon (a square of the
  micro tile's extent, `tests/python/geotiff_fixtures.py`, 3 × 4 cells;
  `@tester` records its byte count), `mesh` with its link as `domain`: exit
  0, one `.vtk` link whose bytes parse with `vtkread`.

## 7. Left out

- **Large uploads over a plain HTTP route** (`PUT /files/<user>/<name>` with
  the bearer token, streamed, larger limit): useful to a client with a shell
  and a big DEM; a no-shell client cannot use it. Question 2.
- **Data URIs inline in any tool argument:** one route (upload, then link)
  is enough; inline content would write files during the check.
- **Another user's link as an input** (sharing): refused. Question 1.
- **Fetching any other URL:** never (Ola's terms).

## 8. Estimate

`python3 tools/count_loc.py`: `workspace.py` about 30, `schema.py` about 12,
`server.py` about 55 (`_inputs` 20, the tool and `_upload` 35), `store.py`
about 15, `http.py` about 10, the README paragraph 0. **About 125; worst
case 180**, one PR, under CLAUDE.md §2's 700. No dependency added: `base64`,
`os` and `hmac` are the standard library; none is on §2's prohibited list.

## 9. Ideas file

Read `../rasputin_scratch/ideas/ideas.md`, 2026-10-19. Takes up the next
step of "The MCP server as a service: results exchanged without a shared
filesystem" (inputs without a shared disk). Overlaps "MCP tool descriptions:
how to mesh a domain too big for one run" only in editing the instructions;
not taken up. Left out and handed back: the large-upload route and link
sharing (section 7).

## 10. Questions for Ola

1. May a user pass a link to another user's file as an input (sharing a
   result between colleagues)? *Default: no; refused in plain words.*
2. Should a later step add a plain HTTP upload for large files, for clients
   that have a shell? *Default: not now; uploads go through the tool, up to
   2 MiB.*

## Review

- **Design review, round 1, 2026-10-19.** On `87db2b97`, 0 counted lines (estimate about 125, at most 180), @reviewer: APPROVED WITH WORDINGS. Safety holds (no address fetched, read values only, relative paths, token-and-user check sound); SDK 4 MiB limit and TiffFile stream reading checked. Five wordings applied by the main session: the link prefix as `links.url` builds it; write-URI refusal HTTP only (and T5); base64 length bound plus T4 cases; two more `fullmatch` killers; a same-length look-alike host.
- **Code review, round 1, 2026-10-20.** On `32ad828a` (red `fcc4401f`, green `8bceeafe`, README `32ad828a`), 142 counted lines (estimate about 125, at most 180), @reviewer: APPROVED. Safety points checked in code; five named mutants killed (red commit message); README and HTTP_INSTRUCTIONS match code; full suite with -rE 8960 passed, 0 errors (the developer's one error did not recur). Suggestions: `_upload` lets a `ValueError` from `Store.receive` (store.py:78, untested) and a text encode error (lone surrogate) reach the client as raw exceptions, the first carrying server paths; the base64 over-limit size overstates by up to 2 bytes.
- **Code review, round 2, 2026-10-20.** On `1fbf6da1` (red `441bd6cc`, green `1fbf6da1`), 149 counted lines (estimate about 125, at most 180), @reviewer: APPROVED. Both upload refusals plain words with no path or exception name; the store's text logged as a warning only; 13 MCP suites pass (328 tests). CI green owed before merge.
