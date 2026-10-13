# Increment 37: an MCP server over the rasputin commands

**Status: built and review-approved** on branch `worktree-mcp-server` (red `ee5b5d52`, `@tester`, 62 tests; green `851592ba`, `@developer`, 322 net production lines by `python3 tools/count_loc.py b4a81c0d HEAD` against the estimate of 260 and the worst case of 360; no mutation round and no `@perf` run, section 10; section 12 is the as-built record; code review round 1 CHANGES REQUESTED, four blockers, all fixed: T5 case and prose `783c3368` (`@tester`), the `error` field `730162fe` (`@developer`, 323 net lines), the README section `19652b88` (`@architect`), two citations pinned `dbc77f3b`; round 2 CHANGES REQUESTED, two prose blockers, fixed: the README wording `2a871c50`, the help-text escape `01025786`; master `c317c7e4` merged in as `030fd32b` (328 net lines by `count_loc.py origin/master HEAD`), the one test it broke fixed `d8a1729d`; round 3 APPROVED; next the push on Ola's yes). Designed by
`@architect`, 2026-10-12, on `b4a81c0d` (branch `worktree-mcp-server`); the
nine defaults of section 6 were taken while Ola was away. Probe and its
output: `docs/increments/37-probes/`. Design review round 1 (Review section)
asked for three fixes, in sections 6.4, 6.5 and 7 and tests T4 to T6, T9 and
T11; after round 2, `@tester`'s red suite (section 9).

## 1. What Ola asked

An external agent (Codex, ChatGPT, Mistral, Claude Desktop) should discover
and call rasputin as tools instead of parsing its help text. Ola, 2026-10-11:
"my working name for what's we're doing is intent-driven triangulation. This
means accessing the tools using agents and intent. [...] How accessible is
the tooling to an agent like mistral or chatgpt?"; "Add an MCP server to the
ideas file after GeoPackage." The dependency is ruled: "Rule the mcp SDK: yes,
use it" (2026-10-12), so the server is built on the official `mcp` Python
SDK.

The thing exposed is the CLI as it is today, one tool per command, with a
JSON schema per tool, structured errors (the plain-word refusal text plus the
exit code), and the run record, `--stats` and output paths returned as data.
The GeoPackage output (increment 36) is the mesh an agent can read. Nothing
in the server reaches the C++ core or widens the I/O boundary: the server
runs the CLI, and the CLI keeps every rule it has.

The commands today (`src_python/tin_engine/cli.py@b4a81c0d:169`, listed by
`typer.main.get_command(app).commands` as `tests/python/cli_driver.py:42`
does): `version`, `draw`, `mesh`, `palette`, `fetch`, `fetch-stations`,
`catchment`, `station-catchments`. Their parameter counts, from the probe:
mesh 34 (36 after increment 26b merged, section 12), draw 10, fetch 10, catchment 11, station-catchments 11,
fetch-stations 3, palette 2, version 0.

## 2. Prior art: legacy and literature

*Literature.* The protocol is the Model Context Protocol, specification
revision **2026-07-28** (modelcontextprotocol.io/specification/2026-07-28;
schema `schema/2026-07-28/schema.ts`), JSON-RPC 2.0 messages over a
transport. What this design uses, by the specification's own sections:

- *Server Features, Tools*: a tool is `name`, `description`, `inputSchema`
  (JSON Schema, 2020-12 by default), optional `outputSchema` and
  `annotations` (`readOnlyHint`, `destructiveHint`, `idempotentHint`,
  `openWorldHint`). A result carries `content` blocks and optional
  `structuredContent`, which must conform to `outputSchema` when one is
  given. Two error channels: *protocol errors* (unknown tool, malformed
  request; JSON-RPC error `-32602`) and *tool execution errors* ("input
  validation errors", "business logic errors"; a result with
  `isError: true`, which "clients SHOULD provide to language models to
  enable self-correction"). Servers "MUST validate all tool inputs". Tool
  names should use letters, digits, `_`, `-` and `.`. A tool "MAY return
  links to Resources" as `resource_link` content.
- *Utilities, Progress*: a client that wants progress sends a
  `progressToken` in `_meta`; the server MAY send `notifications/progress`
  with `progress`, optional `total` and `message`; "the `progress` value
  MUST increase with each notification, even if the total is unknown", and
  it "MAY be floating point"; both sides "SHOULD implement rate limiting".
- *Utilities, Cancellation*: over stdio a client cancels with
  `notifications/cancelled` naming the request id; the server SHOULD stop
  processing, free resources and send no response; both sides handle the
  race where the work finished first. Its *Timeouts* section: a sender
  SHOULD time out a request and cancel it; "implementations MAY choose to
  reset the timeout clock when receiving a progress notification
  corresponding to the request", but "SHOULD always enforce a maximum
  timeout, regardless of progress notifications".
- *Extensions, Tasks* (`io.modelcontextprotocol/tasks`): durable handles for
  long operations, `tasks/get`, `tasks/cancel`, status `working`,
  `completed`, `failed`, `cancelled`. An opt-in extension both sides must
  negotiate; "host support varies by client". The installed SDK's client
  dropped its experimental task support in v2 (migration guide, "Experimental
  Tasks removed") and the low-level server has no task handler, so Tasks are
  out of this increment (section 6.5).

The SDK is `mcp` **2.3.0** (PyPI; MIT; `Requires-Python >=3.10`;
py.typed), measured in the probe venv. Its v2 facts this design rests on
(py.sdk.modelcontextprotocol.io/v2, confirmed by introspection in the probe):
`FastMCP` is renamed `MCPServer`; the low-level `Server` takes handlers as
constructor arguments (`on_list_tools`, `on_call_tool`, each
`async (ctx, params) -> result`) and runs over `mcp.server.stdio.stdio_server()`;
`Client(server)` connects in process for tests, `Client(StdioServerParameters(
command, args, env, cwd))` spawns a stdio server, and `Client("http://...")`
speaks Streamable HTTP; `Client.call_tool(name, arguments,
read_timeout_seconds, progress_callback)`; `ctx.session.send_progress_notification(
progress_token, progress, total, message)`; `ctx.meta["progress_token"]` is
the token (the wire key is `progressToken`); an `MCPError` raised in a
handler becomes a JSON-RPC error, any other exception "Internal server
error"; an async handler's `await` is cancelled when the client cancels
(the SDK's cancellation page), a sync handler is not. Wire fields are
snake_case in Python (`input_schema`, `is_error`, `structured_content`).
Two facts the probe found that the documentation does not say: **the SDK
validates no tool argument on the server side** (its `jsonschema` use is in
the client; `MCPServer` accepted an unknown argument and ran the tool,
probe run 1), and **the in-process transport carries no progress token**
(`mcp/shared/jsonrpc_dispatcher.py` sets it; the in-process
`direct_dispatcher` does not; probe part E, "progress ticks 0" in process,
1 over stdio). The SDK's dependency set is 13 packages (`anyio`, `httpx2`,
`jsonschema`, `mcp-types`, `opentelemetry-api`, `pydantic`, `pyjwt[crypto]`,
`python-multipart`, `sse-starlette`, `starlette`, `typing-extensions`,
`typing-inspection`, `uvicorn`), none on CLAUDE.md §2's list; importing
`mcp.server.lowlevel` and `mcp.client` costs 0.40 s (probe E).

Other geospatial tools over MCP, as found (sparkgeo's tracked list,
`github.com/willcadell/geo-mcp-servers`, "100 servers tracked" on
2026-10-12): `gdal-mcp`
(`github.com/JordanGunn/gdal-mcp`): one tool per operation (`raster_info`,
`raster_reproject`, `vector_clip`, ...), FastMCP 2.0 over stdio, built on
Rasterio, GeoPandas and pyogrio, which this project may not use; the QGIS
plugins (`qgis_mcp_plugin`, 118 tools; `qmcp`) run a socket server inside
QGIS Desktop. None wraps an existing CLI with schemas derived from it; all
hand-write their tool signatures. This design differs in that one respect
only, and claims no novelty: it is an adapter to a published protocol.

Typer's own switch for plain errors, used by the server for the child
process: `TYPER_USE_RICH` (Typer 0.27.3, `core.py` line 27 of the package,
`HAS_RICH = parse_boolean_env_var(os.getenv("TYPER_USE_RICH"), default=True)`).
With `TYPER_USE_RICH=0` a refusal is `Usage: ...`, `Try ... --help`, a blank
line and one `Error: <words>` line (probe; `_TYPER_STANDARD_TRADITIONAL` has
no effect in this version).

*Legacy.* Nothing. The archived tree has no MCP, JSON-RPC or JSON Schema
code:

    $ git grep -n -i -E "mcp|model context protocol|json-?rpc|json schema|jsonschema" legacy-archive -- legacy
    (no files; exit 1)

Its only servers are a three.js web viewer (`legacy/rasputin/web_visualize.py`,
`git grep -l -i -E "server|rpc" legacy-archive -- legacy` also lists
`solar_position.h`, a false hit on "observer"); neither carries across.

*In this tree.* `tests/python/cli_driver.py:42-56` already derives the
command set from the Typer app and strips Rich's boxes from refusals; the
server does the first the same way and avoids the second with
`TYPER_USE_RICH=0`. `tools/bench.py:45-90` runs the CLI as a child through a
`Runner` protocol with an injected `SubprocessRunner`; the server's child
runner follows that shape. `src_python/tin_engine/run_record.py@b4a81c0d:288`
(`as_json`) is the record the server returns.

## 3. Measured on this worktree (`docs/increments/37-probes/mcp_probe.py`)

AC power, Python 3.14.7, mcp 2.3.0, typer 0.27.3, pydantic 2.14.1; the
probe's own venv, `tin_engine` from this worktree with the `_core` built at
`46b5b1cb` (no C++ changed since: `git diff --stat 46b5b1cb HEAD -- include
src bindings CMakeLists.txt` is empty); `@developer` was building in another
worktree, so timings are medians of 3 to 5 on a shared machine. The full
output is `docs/increments/37-probes/mcp_probe_output.txt`.

The probe is the design in miniature: a low-level `Server` with one tool per
Typer command, each schema derived from the command's click parameters,
arguments validated by the derived Pydantic model, the command run as a child
process.

| what | measured |
|---|---|
| schemas derived from all 8 commands | 16,575 bytes of JSON; `mesh` 34 properties, `required: ["out"]`, 8,269 bytes |
| parameters needing a note | 2 of 81 (at `b4a81c0d`; 83 after increment 26b, section 12): `snap_spacing` on `draw` and on `mesh` has a validating callback (`_snap_spacing`, `src_python/tin_engine/cli.py@b4a81c0d:478`, `:676`), which stays in the child and still refuses there |
| tuple options | `bbox` -> `prefixItems` of four numbers, `minItems`/`maxItems` 4; `tolerance_near` -> `[string, number]` |
| flag pairs (`--binary/--ascii`, `--delaunay/--no-delaunay`, `--[no-]features-merge-same-class`) | one boolean each, default from the CLI (`binary: true`) |
| `fetch --connections` (an integer range) | `{"type": "integer", "minimum": 1, "default": 8}` |
| refusal, missing `--out` | caught by the derived model before any child runs: `is_error: true`, `"out: Field required"` |
| refusal from the child (bad `--dem`) | `is_error: true`, `{"exit_code": 2, "error": "exit 2: Invalid value for --dem: cannot read /nonexistent.tif: No such file or directory"}` |
| wrong type, unknown argument | `"name: Input should be a valid string"`; `"bogus: Extra inputs are not permitted"` (`additionalProperties: false`) |
| unknown tool | JSON-RPC error `-32602`, "Unknown tool: nosuch" |
| client timeout 1.0 s on a tool awaiting a child `sleep 30` | client raises after 1.00 s; the child is dead 0.5 s later, in process and over stdio |
| two concurrent `version` calls | 0.29 s together (one alone 0.26 s): they ran side by side |
| tool round trip, `version` (child inside) | in process 0.261 s; over stdio 0.264 s; stdio connect and initialize 0.465 s |
| `rasputin mesh` on the 1 m benchmark (`tests/fixtures/dem_archive/7908_3_10m_z33.tif`, 5051 x 5051 nodes at 10 m, `--tolerance 1`) | child 0.72 s (0.72, 0.72, 0.72); 464,290 triangles; binary `.vtk` 16.9 MB; record 1,583 bytes, 38 keys |
| the same through the tool | 0.715 s in process, 0.719 s over stdio; result carries the record (38 keys), `outputs` with the mesh path, one stderr line; the text block 1,958 characters |
| progress ticks for that run | 0 in process (no token), 1 over stdio (the one stderr line, the summary at the end) |
| where the child's time goes | `python -c pass` 0.019 s; `import tin_engine.cli` 0.247 s; `rasputin version` 0.262 s as a child, 0.0107 s in process (`CliRunner`) |

So a tool call costs about 0.25 s over the command itself, almost all of it
importing `tin_engine.cli` (which imports `_core`) in the child; against a
mesh of seconds to minutes this is noise, against `version` it is the whole
call. Section 6.4 takes the child route anyway and says why.

## 4. Blueprint: components and data flow

```
agent (Codex, Claude Desktop, ...)        stdio: JSON-RPC 2.0 lines
   |  tools/list, tools/call  ---------------------------------->  rasputin mcp
   |                                                                 |
   |                                            tin_engine.mcp.server (L5)
   |                                              - one Tool per Typer command, from
   |                                                tin_engine.mcp.schema (click params ->
   |                                                Pydantic model, extra forbidden, and back to argv)
   |                                              - validates arguments with that model
   |                                              - runs the child, kills it on cancel
   |                                              - shapes the result; maps exit 2 to is_error
   |                                                 |
   |                                                 |  argv, env TYPER_USE_RICH=0
   |                                                 v
   |                                     python -m tin_engine <command> ...   (a child process:
   |                                       cli.py as today; --record to a sidecar          the only
   |                                       file beside --out, read back and removed)
   |                                                 |                                 place _core
   |  <--- result: exit_code, outputs, record, ------+   stdout lines, stderr lines,   runs
          stats, stderr; or is_error + words             exit code, the record JSON
```

Boundaries:

- **The server never imports `_core`, `raster`, `io/` or any pipeline.**
  Its first-party imports are `tin_engine` (for `installed_version`) and
  nothing else; it takes the Typer app as an argument, so the schema layer is
  tested on a toy app without the kernel (section 9, T1 to T4). This is the
  Pybind11 firewall unchanged: `_core` runs in the child, behind the CLI.
- **The CLI is the contract.** The server adds no option, no default, no
  validation of its own beyond the schema the CLI implies; a rule the CLI
  enforces (paths under `--out-parent`, `--domain` needs `--tolerance`, the
  prohibited-dependency and I/O boundaries) holds because the CLI runs. A
  new command or option appears as a tool or a property with no server
  change; the schema and the CLI cannot drift because the schema is read
  from the CLI at start.
- **Data, not execution.** `tin_engine.mcp.schema` is pure: Typer command in,
  Pydantic model and `Tool` out, and arguments in, argv out. No I/O, no
  async. `tin_engine.mcp.server` holds the only I/O (the child, the sidecar
  record) and the only async code.
- **The server's stdout is the wire.** Over stdio the server process writes
  nothing but JSON-RPC to its own stdout: the child's stdout and stderr are
  pipes the server reads, and anything the server has to say for itself
  (an unexpected exception's traceback) goes to its stderr. A stray
  `print` in the server would corrupt the protocol stream, so `server.py`
  has none, and T7 (a real stdio round trip) is where one would show.
- **Async-ready.** The server is `anyio`/`asyncio` throughout: the child is
  `asyncio.create_subprocess_exec`, its streams are read as tasks, a second
  call does not wait for the first unless the jobs cap says so, and the same
  `Server` object could be mounted in Streamable HTTP later without change.

## 5. The tool surface

One tool per command, named as the command with `-` as `_` (`fetch_stations`,
`station_catchments`): the spec allows `-`, but most clients and models
write tool names as identifiers. `rasputin mcp` itself (section 6.8) is
excluded from the list. Tool order is the CLI's (deterministic, as the spec
asks).

Per tool:

- `description`: the first paragraph of the command's docstring, as Typer
  shows it in `rasputin --help`.
- `input_schema`: the derived model's JSON Schema with `title` removed and
  `additionalProperties: false`. Each property carries the option's help
  text as `description`, its default, and its type by this table, which the
  probe exercised on all 81 parameters at `b4a81c0d`: click `str`/`path` -> string;
  `int` -> integer (`IntRange` adds `minimum`/`maximum`); `float` -> number;
  `bool` flags -> boolean (a `--x/--no-x` pair is one boolean; `true` emits
  the first name, `false` the second); a tuple type -> fixed-length array
  (`prefixItems`); `multiple=True` -> array of the item type; a `Choice`
  -> `enum`; an option with default `None` -> nullable, not required;
  required options and arguments -> `required`. Paths are strings: the
  agent gives paths on the machine the server runs on, as the CLI does.
- `output_schema`: one `RunResult` model for every tool (below), so clients
  can validate.
- `annotations`: `version` `read_only_hint=True, idempotent_hint=True`;
  `fetch` and `fetch_stations` `open_world_hint=True` (network); every tool
  `destructive_hint=False` (the CLI never deletes; it refuses to overwrite a
  mesh with a report, `src_python/tin_engine/cli.py@b4a81c0d:1335-1364`).

`RunResult` (Pydantic, `tin_engine.mcp.server`):

| field | meaning |
|---|---|
| `command` | the argv the server ran, as `shlex.join` prints it (what `--record`'s `command` field says, `src_python/tin_engine/cli.py@b4a81c0d:1367`) |
| `exit_code` | the child's exit code (0 here; a refusal is an error result, below) |
| `seconds` | wall time of the child |
| `outputs` | the stdout lines that are paths of files that exist after the run (mesh: the mesh, `--out-edges`, `--stats` file, `--record` file; catchment: the polygon; fetch_stations, station_catchments: directories) |
| `stdout` | every stdout line (fetch's report, palette's preset text, version's two lines) |
| `stderr` | every stderr line (the mesh summary, feature counts, warnings) |
| `record` | for `mesh`: the `--record` JSON as an object (38 keys on the fixture), whether or not the caller asked for a record file; when the caller did not, it is read from a sidecar file beside `out` (section 6.4) and the sidecar is removed |
| `stats` | for `mesh` with `stats` given: the Markdown report text; `stats: "-"` is redirected to a sidecar file beside `out` and returned here rather than mixed into stdout |

The `content` of a successful result is one text block with the same JSON
(the spec's backward-compatible form) plus one `resource_link` per output
(`file://` URI, the file's suffix as `mimeType` where known: `.gpkg`
`application/geopackage+sqlite3`, `.vtk` and `.ply` `application/octet-stream`,
`.geojson` `application/geo+json`, `.json` `application/json`, `.md`
`text/markdown`). A large mesh never travels through the protocol: the
agent reads the GeoPackage with its own tools.

Errors, following the spec's two channels:

| case | result |
|---|---|
| arguments fail the derived model (missing required, wrong type, unknown key) | `is_error: true`; text and `structured_content` `{"exit_code": null, "error": "<field>: <Pydantic's message>; ..."}`; no child runs |
| the child exits 2 (a `BadParameter`, `src_python/tin_engine/cli.py@b4a81c0d:416-1177`) or 1 (`typer.Exit(1)` after a fetch or service error, `src_python/tin_engine/cli.py@b4a81c0d:1527`, `:1566`, `:1569`, `:2560`) | `is_error: true`; `{"exit_code": 2, "error": "exit 2: <the Error: line's words>"}`; `stderr` kept in full |
| the child dies of a signal | `is_error: true`; `{"exit_code": -N, "error": "killed by signal N (<name>)"}` |
| unknown tool | `MCPError(-32602, "Unknown tool: <name>")`, a protocol error, as the spec files it |
| the server's own fault (an exception in the handler) | caught and returned as `is_error: true` with `exit_code: null` and the exception's words; the SDK's "Internal server error" is never the answer |

The refusal words are the child's `Error:` line(s) with `Usage:` and
`Try ... --help` dropped; measured: "Invalid value for --dem: cannot read
/nonexistent.tif: No such file or directory". They are the CLI's words, so
increment 25's plain-output rules carry over for free.

## 6. Decisions, each with its default for Ola

### 6.1 Dependency: an optional extra (default), not required

`pyproject.toml` gains `mcp = ["mcp>=2.3"]` under
`[project.optional-dependencies]` beside `codecs`, `viewer` and `dev`
(`pyproject.toml@b4a81c0d:48-70`). A mesh-only install is unchanged; `pip
install "rasputin[mcp]"` adds the SDK's 13 packages (section 2), two of them
HTTP servers this increment never starts. `tools/check_prohibited_deps.py`
scans every extra's entries (`tools/check_prohibited_deps.py@b4a81c0d:204-216`,
`for extra in project.get("optional-dependencies", {}).values()`), so the
new extra is checked and passes: none of the 14 names is on the list.
`@developer` runs the checker after the edit, and plants `rasterio` in the
new extra once to see it fail (the required-reading rule on verifying a
gate). CI's Python jobs install `.[dev]` today, the `python` matrix job at
`.github/workflows/main.yaml@ce8c364f:261` and the extras job at `:316`, `:337`,
`:351`; all four lines add `mcp` to the extras so the suite runs there;
locally the suite `pytest.importorskip("mcp")`s.
mypy: both `mcp` and `mcp_types` ship `py.typed`, so no override.

### 6.2 Transport: stdio now; HTTP is one later flag, not designed here

Stdio is what Claude Desktop, Codex and Claude Code launch (a command and
arguments in a JSON config, section 6.8). Streamable HTTP is what a hosted
ChatGPT connector needs; it is the same `Server` object behind
`mcp.server.streamable_http_manager`, with authentication and a listening
port to think about, so it is its own small increment if Ola wants ChatGPT
in particular. Nothing here forecloses it.

### 6.3 Schema source: the Typer command's click parameters

Measured against the alternative. The click parameters of all eight commands
give name, type, default, required, help, repetition and tuple arity for 81
of 81 parameters (at `b4a81c0d`; the count moves with the CLI, section 12); the one note is a validating callback, whose refusal still
comes back from the child in words. A Pydantic model per command would be
a second statement of the same 81 parameters (about 34 fields x 3 lines for
`mesh` alone) that drifts the first time someone adds an option to one and
not the other; the drift test it would need is the derivation itself. So:
`tin_engine.mcp.schema.derive(app) -> list[DerivedCommand]`, where a
`DerivedCommand` is the click command, its tool name, its Pydantic model
(`create_model`, `ConfigDict(extra="forbid")`) and `argv(args) -> list[str]`,
which emits only the keys the caller set (`model_fields_set`), so every CLI
default stays the CLI's.

What is lost, by design: option grouping (`rich_help_panel`), `metavar`
(`BOX`, `FILE N`) and the validator callback's rule are not in the schema;
the help text is the whole documentation of a property. The remedy for weak
help text is better help text in `cli.py`, which then serves both readers.

### 6.4 Execution: a child process per call

In-process invocation is 25 times cheaper (0.0107 s against 0.262 s for
`version`) and was rejected for five reasons, each measured or inherent:
(1) a mesh runs C++ for seconds to minutes with the GIL held in stretches;
a sync handler cannot be cancelled (SDK cancellation page), an async one
cannot interrupt a C++ call; the child is killed within 0.5 s of the cancel
(probe C). (2) Two agents' calls, or one agent's two calls, run as two
processes with their own memory and no shared `sys.stdout`; in process the
CLI's `typer.echo` to the process's stdout would have to be redirected per
call, which is process-global. (3) A crash in `_core` (a bounds check, an
allocation failure at basin scale) ends the child and returns an error; in
process it ends the server. (4) Memory from a basin run is returned to the
system when the child exits; the server stays small. (5) The result is the
exact CLI behaviour, including `sys.argv`-derived fields in the record.
The cost, 0.25 s a call, is importing `tin_engine.cli` in the child; it is
paid once per tool call and is below the noise of any mesh.

The child is `[sys.executable, "-m", "tin_engine", *argv]`: the same
interpreter and the same installed `tin_engine` as the server, found without
`PATH`. That needs `src_python/tin_engine/__main__.py` (new, 3 lines:
`sys.argv[0] = "rasputin"` so the record's `command` field reads
`rasputin mesh ...`, then `app()`); `python -m tin_engine version` becomes a
second way to run the CLI. The child's environment is the server's plus
`TYPER_USE_RICH=0`, `NO_COLOR=1`, `TERM=dumb`; `RASPUTIN_DATA` and the cache
default pass through unchanged.
The child is started without `cwd=`, so it inherits the server's working
directory: a relative `out` and its sidecar beside it resolve to the same
place in the child and in the server's read-back and `unlink`.

*The sidecar files.* For `mesh` without `record_path`, the server appends
`--record <sidecar>` and reads it back; with `stats: "-"` it appends
`--stats <sidecar>` instead and returns the text; both sidecars are removed
after the read (`unlink(missing_ok=True)`, in the handler's `finally`, so a
refused or killed run leaves none). The sidecar is **beside `out`**, not in
the system temp directory: `Path(out).parent / f".{Path(out).name}.record-
{token}.json"` and `...stats-{token}.md`, `token` eight hex digits from
`secrets.token_hex(4)`, the name chosen by the server and the file created
by the child. Checked against `_destination`
(`src_python/tin_engine/cli.py@b4a81c0d:420-447`), which every `--record`
and `--stats` path passes through (`_report_target` at `:1335-1347`,
`_record_target` at `:1350-1364`, both called at `:1195-1196` after `--out`
at `:1180`):

- `out.is_symlink()` (`:438`): a fresh random name is not a symlink.
- containment in `out_parent` (`:441-444`): `Path(record).resolve()` is
  inside `out_parent.resolve()` exactly when `Path(out).resolve()` is, since
  the two share a directory and `out` is not a symlink (refused at `:438`
  before this); so a system temp path (`/var/folders/...`) would be refused
  whenever `out_parent` is set, and the sidecar beside `out` never is.
- `target.parent.is_dir()` (`:445`): the sidecar's directory is `out`'s,
  which `_destination(out, ...)` at `:1180` has already accepted, or refused
  first naming `--out`; the server never creates the directory or the file.
- `target in taken` (`:1343`, `:1360`): the random name collides with no
  mesh, edges or report path.

The sidecar is created by the child with the child's permissions in a
directory the caller already chose to write a mesh into, so it widens no
boundary. A caller who wants the record kept passes `record_path`, which
is then the caller's path and not removed. The runner is injected (`Child`
protocol: argv prefix in, a `(returncode, stdout lines, stderr lines)`
out), so tests substitute a fake child (section 9).

### 6.5 Long runs: one blocking call with progress; cancellation kills the child; no Tasks

A basin mesh takes minutes and the call stays open for them. Over stdio the
server sends `notifications/progress` when the client gave a token: one per
stderr line the child prints (message: the line), and a heartbeat every 5 s
with message `"running, N s"` when no line came, so a client that resets
its timeout on progress (the spec's Cancellation page allows it) does not
give up on a silent run. Measured: `mesh` prints nothing until its summary
line (probe E: one stderr line in 0.72 s), so on today's CLI the heartbeat
is the only progress a long mesh shows; feature and tolerance runs print a
few lines as phases pass (`src_python/tin_engine/cli.py@b4a81c0d:1084`,
`:1708`, `:1783`, `:1810`, `:1916`).

*The `progress` value.* The spec says it MUST increase with each
notification; a line count would not, since a heartbeat has no line. So
every notification of a call, line or heartbeat or waiting (section 6.6),
carries `progress` = seconds elapsed since the handler began, from
`time.monotonic()`, rounded to milliseconds, and `total` omitted (the run's
length is unknown). One sender per call holds the last value sent and sends
`max(elapsed, last + 0.001)`, so two lines read in the same millisecond
still increase; the 0.001 is the rounding step, and the sender is one
`asyncio` task per call (the stderr reader and the heartbeat timer feed the
same sender), so the values are never sent out of order. The probe's
version (`docs/increments/37-probes/mcp_probe.py:186-188`) sent the line
count; its measured output, the count and text of the ticks, is unchanged
by this choice, so the probe was not rerun.

The constant 5 s assumes runs from seconds to an hour (12 to 720
notifications) and was checked at the fixture's 0.72 s run (zero
heartbeats) and the probe's 1 s cancellation; it is a module constant, not an
option. A client's own timeout is the client's; the SDK's client default is
none.

Cancellation: the SDK cancels the handler's `await` on `notifications/cancelled`
(stdio) or on a dropped stream; the handler's `finally` kills the child and
waits for it, then returns nothing (the spec says no response).

*What a killed run leaves on disk.* The CLI writes a GeoPackage as
`<target>.part` and renames it at the end
(`_write_geopackage`, `src_python/tin_engine/cli.py@b4a81c0d:1250-1289`,
the `os.replace` at `:1289`),
so a killed `.gpkg` run leaves at most a `.part` file, which the next run
replaces. `.vtk` and `.ply` are not written that way: the file is encoded in
memory and written with one `target.write_bytes(data)`
(`src_python/tin_engine/cli.py@b4a81c0d:1241`), and `_writing` (`:2121`)
only turns an `OSError` into a refusal, it does not rename. So **a `mesh`
call cancelled during that write can leave a truncated `.vtk` or `.ply`
under the asked name**; cancelled before it, the old file (if any) is
untouched, and cancelled after it, the mesh is complete and only the report
or record is missing. The server does not guard against this: it cannot
tell a truncated file from a complete one or from the caller's older file
under that name, and deleting on cancel would destroy the last two. It
accepts and documents it instead: a cancelled call returns no result, so
the agent holds no path it was told to trust, and the next run overwrites
the file. The remedy is in `cli.py`, every output through `.part` and a
rename; it is in Ola's ideas file and is not this increment's. The server's
`instructions` (section 6.8) gain one line: "a cancelled `mesh` may leave an
unreadable `.vtk`/`.ply` at `out`; rerun before reading it".

The Tasks extension is the right shape for hour-long runs with
reconnects, and is out: the SDK's v2 client removed its experimental task
API, the low-level server has no task handler, and client support "varies".
When the SDK ships it, a task is a wrapper round the same child runner.

### 6.6 Concurrency: a jobs cap, default 1

Two calls run as two children (probe D), and two basin runs are two to six
GB of memory (`docs/increments/23-basin-scale.md`'s limits). `rasputin mcp
--jobs N` (default 1) is an `anyio.Semaphore`; a call beyond the cap waits,
with heartbeats every 5 s saying `"waiting for a running job, N s"` and the
same elapsed-seconds `progress` as section 6.5 (the clock starts when the
handler begins, before the wait, so the first running notification is above
the last waiting one), and is still cancellable while waiting. `version` and `palette` also wait: a special
case for cheap commands is a second rule to explain for a 0.3 s gain. The
default was checked at two `version` calls; an agent wanting parallel
fetches raises it.

### 6.7 Resources: none listed; outputs are resource links and paths

No `resources/list`: a server that lists every mesh it ever wrote needs a
registry and a retention rule, and a 17 MB to 1 GB mesh does not belong in a
context window. The run record is in the result; each output file is a
`resource_link` with a `file://` URI and a path in `outputs`, which the agent
reads with its own tools (sqlite3 on the GeoPackage, the Gemini note's
`rasputin://meshes/{id}` is the idea this does not take up).

### 6.8 Where it lives and how it starts

`src_python/tin_engine/mcp/` with `__init__.py` (empty), `schema.py` and
`server.py`; `tests/python/test_layering.py` gains rows in L5 (`cli` row adds
`mcp.server`; `mcp` ``; `mcp.schema` ``; `mcp.server` `mcp.schema tin_engine`;
`__main__` `cli`), all sideways or down; `test_the_table_names_every_module`
(`tests/python/test_layering.py@b4a81c0d:143-149`) demands them.

Start: `rasputin mcp [--jobs N]`, a Typer command in `cli.py` that imports
`tin_engine.mcp.server` inside the function and, when the import fails,
refuses with exit 2 and `Error: the mcp extra is not installed; run
pip install "rasputin[mcp]"`. A subcommand is chosen over a console script
because it is where an agent reading `rasputin --help` finds it, and the
install story stays one package. The server excludes its own command from
the tool list. A client configuration, which README gains in a short
"Connect an agent" section with the Claude Desktop and Codex forms:

```json
{"mcpServers": {"rasputin": {"command": "/path/to/venv/bin/rasputin", "args": ["mcp"]}}}
```

The server's `instructions` (the initialize result's field, shown to the
model) are six lines: paths are on this machine; `mesh` writes a file and
returns its record; write `.gpkg` to inspect a mesh; `fetch` first for a
catalogue source; a refusal's words say what to change; a cancelled `mesh`
may leave an unreadable `.vtk`/`.ply` at `out` (section 6.5).

### 6.9 Tool names: underscores

`fetch-stations` is tool `fetch_stations`. Default yes; the alternative,
the hyphenated name, is spec-legal but some clients reject it.

## 7. Degeneracy policy

| input | behaviour |
|---|---|
| malformed arguments (missing required, wrong type, unknown key, wrong arity for `bbox`) | refused by the derived model before any child runs; `is_error` with Pydantic's words per field (measured, section 3) |
| a missing data folder (`--cache` default with `RASPUTIN_DATA` unset, a `--dem` directory that is not there) | the child refuses with its words and exit 2; the server passes them through; it never guesses a path |
| a refused command (any `BadParameter`, any `typer.Exit(1)`) | `is_error: true`, exit code and words; stderr in full |
| a killed run (client cancel, client gone, server stopped) | the child is killed and waited for within the cancel; no response; the sidecar is removed; a `.gpkg` run leaves at most a `.part` file, a `.vtk`/`.ply` run killed during its one `write_bytes` leaves a truncated file under the asked name (section 6.5) |
| `mesh` with `out_parent` set | the sidecar record or report is beside `out`, so `_destination` accepts it exactly when it accepts `out`; the call succeeds, or is refused naming `--out`, never naming `--record` or `--stats` |
| a child that prints a non-UTF-8 byte | decoded with `errors="replace"`; never a server fault |
| a child whose stdout path line is not a file (fetch's report lines) | kept in `stdout`, not in `outputs` |
| two concurrent calls | both accepted; the second waits for the semaphore (jobs 1) with heartbeats, or runs beside the first (jobs 2+) |
| the client never sends a progress token | no notifications; the heartbeat task is not started |
| a tool call for `mcp` itself | not a tool: unknown tool, protocol error |
| the `mcp` extra absent | `rasputin mcp` exits 2 with the install line; every other command unchanged |
| an exception inside the server's own code | `is_error: true` with the words, `exit_code: null`; the server keeps serving |

## 8. Left out, with reasons

- **An intent parser (`parse_mesh_intent`) and a `MeshIntent` schema.** The
  agent is the intent parser; rasputin gives it the typed tools. The Gemini
  note's version (unverified, in Ola's private ideas folder) is a second
  schema beside the CLI's; this design has one.
- **An OpenAPI export.** Cheap once `schema.py` exists (the same models'
  JSON Schema under an OpenAPI 3.1 envelope), and useful for ChatGPT Actions
  and function calling; its own small increment, after this one proves the
  derivation. Handed back as an idea.
- **Streamable HTTP, authentication.** Section 6.2.
- **Tasks.** Section 6.5.
- **Listed resources, a mesh registry.** Section 6.7.
- **Phase progress from `mesh`.** The CLI is silent until its summary; a
  `--progress` or verbose phase line per `PhaseClock` phase would make the
  heartbeat informative. Handed back as an idea; a `cli.py` change for
  humans and agents alike, not an MCP change.
- **Rewriting tool descriptions.** They are the commands' docstrings, which
  reference increments ("increment 16c"); the product review in Ola's ideas
  file covers that rewrite, and it will reach the tools through the
  derivation.

## 9. Tests (`@tester`), red first

All in `tests/python/test_mcp_*.py`, `pytest.importorskip("mcp")` at module
top; the SDK's `Client` runs the server in process except where stdio is
the point. Fixtures and their figures as measured in section 3.

- **T1 schema from a toy app** (no kernel): a Typer app with two commands
  covering every row of section 5's type table (str, path, int, IntRange,
  float, flag pair, 4-tuple, 2-tuple of path and float, `multiple`, Choice,
  nullable, required argument) derives the expected properties, `required`,
  `additionalProperties: false`, defaults and descriptions; the tool name
  maps `-` to `_`; the server's own command is excluded.
- **T2 argv round trip on the toy app**: for each case, `argv(model(args))`
  is the expected command line: only set keys appear; a `false` flag emits
  the second name; a tuple expands; a list repeats the option; an argument
  is positional.
- **T3 the real app**: `derive(app)` yields the eight commands with their
  parameter counts (34, 10, 10, 11, 11, 3, 2, 0; as built the suite pins the
  tool names and the equality with the click parameters instead, section 12)
  and `mesh`'s `bbox`,
  `tolerance_near`, `binary`, `features`, `snap_spacing` schemas as section
  3 shows; the whole set under 20,000 bytes (16,575 measured) so an agent's
  tool list stays small.
- **T4 validation before the child**: missing `out`, `tolerance: "ten"`, an
  unknown key, `bbox` of three numbers: `is_error`, the field named, and the
  fake child never called. And the argv the fake child receives for `mesh`
  with `out: <dir>/a.vtk` and `out_parent: <dir>`: `--out-parent <dir>` is
  passed through, `--record` names a file in `<dir>` whose name starts with
  `.a.vtk.record-` and ends `.json`; with `stats: "-"` the `--stats` value
  is `<dir>/.a.vtk.stats-<8 hex>.md` and no `--record` is added when
  `record_path` is given.
- **T5 refusal mapping**: a fake child returning exit 2 with
  `Usage:`/`Try`/`Error: Invalid value for --dem: ...` lines gives
  `{"exit_code": 2, "error": "exit 2: Invalid value for --dem: ..."}`; exit
  1 likewise; a signal exit gives `killed by signal`. After a fake child that
  exits 2 without creating the sidecar, and after one that creates it and
  exits 2, no `.a.vtk.record-*` file remains in `<dir>`.
- **T6 a real `mesh`** through the tool on the micro tile
  (`tests/python/geotiff_fixtures.py:97-105`, 3 x 4 cells at 10 m x 5 m,
  no codec needed): exit 0, `outputs` is the `.vtk`, `record` is a dict with
  the record's keys (`run_record.WORDING`'s order), the sidecar record file
  is gone (the directory holds the `.vtk` and nothing else), the `.vtk`
  parses with `vtkread`; the same call with `out_parent` set to the `.vtk`'s
  directory succeeds the same way (the case review round 1 found: a sidecar
  in the system temp directory would exit 2 here with "is outside the
  permitted parent"), and with `out_parent` set to a sibling directory it is
  `is_error` naming `--out`, not `--record`; with a relative `out` under
  `monkeypatch.chdir` the mesh and the sidecar's removal land in that
  directory (the child inherits the server's cwd); `outputs` is exactly the
  `.vtk`, so the echoed sidecar path (`src_python/tin_engine/cli.py@b4a81c0d:1247`) is excluded;
  with `stats: "-"` the result's
  `stats` starts with the report's heading and stdout carries no Markdown;
  and once on the Kartverket fixture under `needs_codecs` with
  `--tolerance 1`: 464,290 triangles in the record, 0.7 s on the probe's
  machine.
- **T7 stdio round trip**: `Client(StdioServerParameters(sys.executable,
  ["-m", "tin_engine", "mcp"]))` lists the eight tools and calls `version`;
  the result's `stdout` has two lines and `exit_code` 0. Connect and
  initialize measured at 0.465 s, so the test's budget is generous.
- **T8 cancellation**: the fake child is `sys.executable -c "import time;
  time.sleep(30)"` with its pid written to a file; `read_timeout_seconds=0.5`
  raises at the client and the pid is gone within 1 s (measured: dead 0.5 s
  after a 1.0 s timeout, both transports); over stdio as well, since that is
  the transport the cancel notification travels on.
- **T9 progress over stdio**: a fake child printing three stderr lines 0.2 s
  apart gives three progress notifications with the lines as messages; the
  heartbeat constant monkeypatched to 0.1 s on a child that prints one line
  at 0.25 s and exits at 0.6 s gives at least three `"running, N s"`
  heartbeats around that line; over the whole sequence, lines and heartbeats
  together in the order received, every `progress` is strictly greater than
  the one before, each is within 0.05 s of the wall time since the call
  began (so the value is elapsed seconds, not a count), and `total` is
  absent. With jobs 1 and a first call holding the semaphore, the second
  call's `"waiting for a running job, N s"` heartbeats and its later running
  notifications form one strictly increasing sequence. In process no token
  arrives (the SDK's gap, measured) so the in-process case asserts only that
  the run completes.
- **T10 jobs cap**: with jobs 1, two calls on a 0.5 s fake child complete
  in about 1 s with non-overlapping child intervals; with jobs 2 in about
  0.5 s.
- **T11 `rasputin mcp` without the extra**: the import made to fail
  (`monkeypatch.setitem(sys.modules, "mcp", None)`, and
  `monkeypatch.delitem(sys.modules, "tin_engine.mcp.server", raising=False)`
  so an earlier test's import of the server does not satisfy the lazy
  import), exit 2 and the install line; `rasputin --help` still lists
  `mcp`.
- **T12 layering**: the new rows (section 6.8) pass `test_layering.py`
  unchanged in rule; `test_the_table_names_every_module` sees `__main__`,
  `mcp`, `mcp.schema`, `mcp.server`.

Mutation testing is not required: no suite here is invariant-critical in the
README's sense (no topology); the one place a typo hides, the type table of
6.3, is pinned by T1's exhaustive cases.

## 10. Estimate, speed, acceptance

| file | net lines (counted: Python, not tests or docs) |
|---|---|
| `src_python/tin_engine/mcp/schema.py` (type table, `DerivedCommand`, model, argv, `Tool`) | about 85 |
| `src_python/tin_engine/mcp/server.py` (`RunResult`, child runner with progress, heartbeat and kill, the sidecar names and their removal, the increasing-progress sender, result shaping, errors, semaphore, `serve()`) | about 150 |
| `src_python/tin_engine/mcp/__init__.py`, `src_python/tin_engine/__main__.py` | 3 |
| `src_python/tin_engine/cli.py`: the `mcp` command with `--jobs`, the lazy import and its refusal | about 20 |
| `pyproject.toml`, `.github/workflows/main.yaml`, `README.md`, `test_layering.py` rows | 0 (not counted) |
| total | about 260 (250 before review round 1; the sidecar naming and the progress sender add about 10); with the +39 % worst case of increment 33, about 360: one PR, under 700 |

Speed of existing commands: expected unchanged. The increment adds a
command to the Typer app (registration at import, microseconds), a
`__main__.py` nobody else runs, and a new package nothing else imports; no
refine or mesh code is touched, so no `@perf` acceptance run is due. The
probe's `rasputin mesh --tolerance 1` on the benchmark tile, 0.72 s, is the
same CLI the server will run.

## 11. Ideas file

Taken up: "MCP server over the rasputin commands" (this increment, as ruled).
Overlapped: "Intent layer over rasputin" (its MCP server step is this one,
built from the CLI rather than from a `MeshIntent` schema; its OpenAPI
export and provider benchmark are not taken up); "Dual-use applications" and
"Feature-adaptive vertical tolerance" (the direction, intent in and terrain
out, is what the tools are for; nothing of them is built here); "Release
hygiene with priority" (README's "Connect an agent" section lands here, the
rest stays there); "Product review: readable code docs" (tool descriptions
inherit the docstrings it will rewrite). Handed back under *Ideas* in the
handback: phase progress lines from `mesh`, an OpenAPI export from
`schema.py`, Streamable HTTP for hosted agents, Tasks when the SDK has them,
a jobs cap by memory estimate rather than count.

## 12. As built (recorded by the main session from the handbacks)

Ola picked the MCP server on 2026-10-11 ("Add an MCP server to the ideas
file after GeoPackage") and ruled the dependency on 2026-10-12 ("Rule the
mcp SDK: yes, use it"); the design's nine questions were decided on their
defaults while he was away (section 10). Design `b990e38c`, three review
rounds (Review section), the last two closed on the reviewer's exact
wording in the recording commit.

**Red** `ee5b5d52` (`@tester`): 62 tests, every one red on the missing
module `tin_engine.mcp`, the missing `tin_engine.__main__`, the stdio
server dying at import, or `rasputin mcp` not being a subcommand; the
existing 7,233 tests pass (full suite 9 min 56 s, coverage 98.6 %).
Files: `tests/python/test_mcp_schema.py` (T1-T3, 26),
`test_mcp_server.py` (T4-T6, T8, T10 in process, 22), `test_mcp_stdio.py`
(T7-T9, 5), `test_mcp_cli.py` (T11, 2, run without the SDK too),
`test_layering.py` rows (T12, 7), helpers `mcp_driver.py`,
`mcp_fake_child.py`, `mcp_stdio_launcher.py`. The three SDK suites skip
when `mcp` is absent. Pins beyond the design, asserted: `derive(app)` drops
the command named `mcp`; `DerivedCommand.tool_name`, `.model`, `.argv()`,
`.tool()`; `argv[0]` is the hyphenated CLI name and argv follows the
command's parameter order; an Enum default is its value, a Path default its
string; `build_server(app, *, jobs=1, child=<argv prefix>)`; the
`HEARTBEAT_SECONDS` module constant; error results carry `stderr` in
`structured_content`; the success text block is the structured JSON
verbatim. Departures from the design: `__main__.py` carries an
`if __name__ == "__main__"` guard (the layering test imports it); the
Kartverket triangle count is a window (410k to 520k) with the tolerance read
from the record, since the count is not a record field; "within 0.05 s" on
progress became exact two-sided bounds.

**Green** `851592ba` (`@developer`): `src_python/tin_engine/mcp/`
(`__init__.py`, `schema.py` 98 lines, `server.py` 204), `__main__.py` (5),
`cli.py` (15: the `mcp` subcommand and the `_destination` hint),
`pyproject.toml` (the `mcp` extra, `mcp>=2.3`, pulled by `dev`). No test
file touched. All 62 red tests green; full `pytest tests/python/` 7,292
passed, 17 skipped, coverage 98.43 % (`mcp/schema.py` 100 %,
`mcp/server.py` 91 %); Kartverket through the tool 464,290 triangles in
0.77 s; `ruff check`, `ruff format --check`, `mypy` (strict, 72 files),
`check_prohibited_deps.py` (fails naming rasterio when planted in the `mcp`
extra, restored) and `check_detria_boundary.py` clean. No C++ touched.
Departures, for the code review to judge: (1) `cli.py`'s `_destination`
gains a `hint` (default `--out`) used as `param_hint` on its three
refusals, and the `--out-edges`, `--stats`, `--record` and `--out-dir`
callers pass theirs: section 6.4 and T6 say an `out` outside `out_parent`
is "refused naming `--out`", but at `b4a81c0d` the CLI said
`Invalid value: X is outside the permitted parent Y`, naming no option;
refusals now read `Invalid value for --out: ...` (Ola's question: keep?
default keep, it is what the design and the plain-output rule ask). (2)
`RunResult` lives in `schema.py`, not `server.py` as section 5 sites it:
`tool()` needs it and the layering rows forbid `schema` importing
`server`; `server` re-exports it. (3) Section 6.1's four `main.yaml`
install-line edits are not made: `tests/python/test_ci_changes.py` pins
those strings; `dev` pulls `mcp>=2.3` instead (the main session's default).

**Size:** 322 net production lines (`python3 tools/count_loc.py b4a81c0d
HEAD`): schema 98, server 204, cli 15, `__main__` 5; the design estimated
260 with a worst case of 360. Under 700.

**README** (`@architect`, code review round 1, blocker 1): section 6.8's
"Connect an agent" section landed in `README.md` between "Output and
viewing" and "How it is built" (`git log --oneline -1 -- README.md` names
the commit), human reader first as Ola ruled on 2026-10-11, about 240 words
of prose plus the Claude Desktop and Codex configuration blocks; every claim
in it was checked against `server.py`, `schema.py`, the `mcp` subcommand and
`pyproject.toml` at `730162fe`, the derived tool list in the worktree venv
(eight tools, `fetch_stations` and `station_catchments` with underscores),
and the two clients' own configuration pages.

**Decided while Ola was away (defaults):** the nine design questions
(section 10); the `dev` extra pulls `mcp` so the `python` CI job runs the
suite; the `_destination` hint kept pending Ola's word.

**After the master merge** (`030fd32b`, @developer): master `c317c7e4`
brought increment 26b and the ruling record; `ROADMAP.md` (row 36 against
row 37) and `tests/python/test_layering.py` (26b's rows against this
branch's, one shared `cli` row importing both `land_cover_fractions` and
`mcp.server`) were resolved by hand; `_core` rebuilt in `build-wt/`; all
gates green, ctest 1031/1031, 328 net lines (the five are `01025786`'s
help-text escape, round 2). One test broke: `test_mcp_schema.py` pinned
`mesh` at 34 parameters, now 36 with `--land-cover` and
`--land-cover-year` (83 across the eight tools; the design's 81 held at
`b4a81c0d`). `@tester` (`d8a1729d`) dropped the per-tool counts rather than
update them: the same test's equality check, schema properties against the
command's click parameters, implies the count, and the pin only detected
options added elsewhere; the eight tool names stay pinned. This changes the
design's T3 wording (section 9) from counts to names plus equality; a probe
dropping `land_cover` from the derived `mesh` schema turns the surviving
check red. The help-text escape (`01025786`): one constant `INSTALL_MCP`,
escaped to `\[mcp]` only when `app.rich_markup_mode` is set, so both
renderings read `pip install "rasputin[mcp]"`.

**Lessons for the next retrospective:** a design claim about existing CLI
wording is checked by running the command before a test pins it (T6's
"naming `--out`"); no worktree had a current `_core` build, so `@tester`
copied one from a merged worktree with identical C++ rather than build
beside `@perf`. A pinned snapshot of another module's shape (a CLI's option
count) in a derivation suite is a resolved value where a rule belongs; pin
the rule and the identity set, never the count.

## Review

Design review round 1, 2026-10-12, @reviewer (`b4a81c0d..b990e38c`, 0 production lines, estimate 250 with worst case 350): CHANGES REQUESTED, three blockers. (1) /Users/skavhaug/projects/rasputin/.claude/worktrees/mcp-server/docs/increments/37-mcp-server.md@b990e38c:382-386 and :450 say every output goes through `.part` and a rename; that holds for the GeoPackage only, `.vtk` and `.ply` are written with `target.write_bytes` (`src_python/tin_engine/cli.py@b4a81c0d:1241`) and `_writing` (:2121) wraps OSError, it does not rename; restate that a killed `.vtk`/`.ply` run can leave a truncated file under the asked name. (2) /Users/skavhaug/projects/rasputin/.claude/worktrees/mcp-server/docs/increments/37-mcp-server.md@b990e38c:354-357 and :255: with `out_parent` set, the server's temporary `--record`/`--stats` path is refused by `_destination` (`src_python/tin_engine/cli.py@b4a81c0d:442-445` via :1360), so every such `mesh` call exits 2; say where the temporary file goes (inside `out_parent`, or beside `out`) and add the case to T4-T6. (3) /Users/skavhaug/projects/rasputin/.claude/worktrees/mcp-server/docs/increments/37-mcp-server.md@b990e38c:364-367 and :396-398: the heartbeat and "waiting" notifications carry no `progress` value, which the specification says MUST increase per notification; define it (elapsed seconds) and have T9 assert the increase across lines and heartbeats. Prose fixes for the same pass: :148 "1 of 81" is 2 (`draw` and `mesh` both have `snap_spacing`); :54 the timeout-reset quote is on the specification's Cancellation page, not Progress; :296 the `python` CI job (`.github/workflows/main.yaml@ce8c364f:261`) also installs `.[dev]`; :97 the server list now says 100. Suggestions: state that the server process writes nothing but JSON-RPC to its own stdout; T11 also drops `tin_engine.mcp.server` from `sys.modules`. Verified: every probe figure against the saved output; SDK 2.3.0 facts in a fresh venv (low-level `Server`, `jsonschema` only in the client session, no progress token on the direct dispatcher, no server-side Tasks, 13 non-Windows dependencies); the specification quotes, the legacy grep, check_citations, and the prohibited-deps scan of extras. Ideas for Ola's file: write every mesh output through `.part` plus rename in `cli.py`; `destructive_hint=False` while `--out` overwrites an existing mesh deserves a second look.

Design review round 2, 2026-10-12, @reviewer (`b9a692de..cb534b36`, 0 production lines): CHANGES REQUESTED, one blocker, applied in this recording commit as the reviewer specified: /Users/skavhaug/projects/rasputin/.claude/worktrees/mcp-server/docs/increments/37-mcp-server.md@cb534b36:363-370 started the child without saying it inherits the server's working directory (no `cwd=`), which the read-back and `unlink` of a relative `out`'s sidecar depend on, and no test pinned a relative `out`; one sentence added there, and T6 gains a relative-`out` case under `monkeypatch.chdir`. The reviewer's idea applied with it: T6 says `outputs` is exactly the `.vtk`, so the echoed sidecar path (`src_python/tin_engine/cli.py@b4a81c0d:1247`) is excluded. Verified: round 1's fixes against `cli.py@b4a81c0d` (`write_bytes` :1241, `_write_geopackage` :1250-1289 with `os.replace` :1289, `_destination` :438-445, the sidecar accepted exactly when `out` is, absolute, relative or through a symlinked parent); the specification quotes on the 2026-07-28 Progress and Cancellation pages; `main.yaml` installs `.[dev]` in four jobs; check_citations exit 0. Lesson: an assumption named only in a handback is not in the design; the file must carry it.

Design review round 3, 2026-10-12, @reviewer (`cb534b36..37f811c1`, 0 production lines): CHANGES REQUESTED, one blocker, applied in this recording commit as specified: /Users/skavhaug/projects/rasputin/.claude/worktrees/mcp-server/docs/increments/37-mcp-server.md@37f811c1:612 cited `cli.py` line 1247 at `b4a81c0d` without its directory path, so `tools/check_citations.py` exited 1 on the branch; now `src_python/tin_engine/cli.py@b4a81c0d:1247`, the record-path echo, so the claim holds. The cwd sentence and both T6 insertions read in place; round 2's citation resolves. The design proceeds to the red suite on this path fix without a fourth round, as the recording-commit practice allows (retrospective 2026-10-12, P1). Lesson for the main session: run `python3 tools/check_citations.py` from inside the worktree after each recording commit; the main checkout's copy judges the main checkout's tree.

**Code review, round 1, 2026-10-12.** Range `b4a81c0d..3859b138`, 322 net production lines (schema 98, server 204, cli 15, `__main__` 5; the overrun over 260 is `server.py`'s shielded clean-up, heartbeat and resource links; under 700), @reviewer: CHANGES REQUESTED, four blockers. (1) /Users/skavhaug/projects/rasputin/.claude/worktrees/mcp-server/docs/increments/37-mcp-server.md@3859b138:508-509 and :665 say the README gains a "Connect an agent" section; `README.md` is unchanged and section 12 records no departure; the section is written by `@architect`. (2) /Users/skavhaug/projects/rasputin/.claude/worktrees/mcp-server/src_python/tin_engine/mcp/server.py@3859b138:147-148 joins every non-Usage stderr line into `error`, where the design (:291-292) says the `Error:` line's words; a `fetch` failing mid-run prints progress lines first (`src_python/tin_engine/cli.py@3859b138:1543`, then :1549), so `error` reads "exit 1: 1,000 of 2,000 bytes Error: ..."; fix: keep the `Error:` lines when any exist, else all; `@tester` adds the T5 case first, `@developer` fixes. (3) Two citations this branch shifted, pinned at `b4a81c0d` in this recording commit: /Users/skavhaug/projects/rasputin/.claude/worktrees/mcp-server/docs/increments/37-mcp-server.md@3859b138:500 (`test_layering.py` 143-149, the function now at 146-152) and /Users/skavhaug/projects/rasputin/.claude/worktrees/mcp-server/docs/research/generic-harness.md@3859b138:192 (`pyproject.toml` 108, now 115). (4) Red-step prose remains in /Users/skavhaug/projects/rasputin/.claude/worktrees/mcp-server/tests/python/mcp_driver.py@3859b138:10-12 and /Users/skavhaug/projects/rasputin/.claude/worktrees/mcp-server/tests/python/mcp_stdio_launcher.py@3859b138:33-34 ("before the increment is built"); `@tester` rewords to the present. Verified: 159 tests green in the worktree venv (Kartverket 464,290 triangles, 0.83 s); sections 5 to 7 against the code (schema from click parameters with `extra="forbid"`, `mcp` dropped; child argv, three environment variables, no `cwd=`; sidecars beside `out`, unlinked shielded; progress strictly increasing with the heartbeat; cancel kills; the semaphore; `file://` links; no print); the I/O boundary (the server imports only `installed_version` and `schema`); prohibited-deps and citation gates exit 0; the green commit touched no test; the three departures sound (each `_destination` hint matches its caller's option; `RunResult` re-exported; `dev` pulling `mcp` reaches every `.[dev]` CI job). Question for Ola: keep the `_destination` hint, default keep. Ideas (for Ola's file): `schema._annotation` drops a `min=0`/`max=0` bound (`if v`); await the cancelled heartbeat task so no tick races the result. Lesson: a design that names a README section needs the as-built record to say it landed or why not.

**Code review, round 2, 2026-10-12.** Range `3859b138..dad373fe`, 323 net production lines, @reviewer: CHANGES REQUESTED, two blockers, both prose. (1) /Users/skavhaug/projects/rasputin/.claude/worktrees/mcp-server/src_python/tin_engine/cli.py@dad373fe:278: the help text `pip install "rasputin[mcp]"` renders as `pip install "rasputin"` under Typer's default Rich markup, which eats `[mcp]` as a tag (measured with `python -m tin_engine mcp --help`; correct only with `TYPER_USE_RICH=0`); `@developer` escapes the bracket for Rich and re-renders, in the same step as the master merge. (2) /Users/skavhaug/projects/rasputin/.claude/worktrees/mcp-server/README.md@dad373fe:135-136 said the server reports every line the command prints; only stderr lines become progress (`src_python/tin_engine/mcp/server.py@dad373fe:116-119`), stdout returns in the result; fixed in this recording commit to "every line the command prints to stderr". Verified: round 1's four fixes against their sources; the two T5 cases match `_refusal` as fixed; 57 MCP tests pass (Kartverket 464,290 triangles, 0.83 s); every README claim against the code and both clients' pages; no red-step scaffolding; check_citations exit 0. Idea (for Ola's file): `python -m tin_engine` usage lines read `python -m tin_engine.rasputin` from the argv[0] rename; `app(prog_name="rasputin")` would fix it. Lesson: help text is a rendered claim; check it as the default renderer shows it.

**Code review, round 3, 2026-10-12.** Range `dad373fe..0615aff7`, 328 net production lines by `count_loc.py origin/master HEAD` (323 at `dad373fe`; the five are the help-text escape), @reviewer: APPROVED, no blockers. Verified: both help renderings print `pip install "rasputin[mcp]"` (Rich default and `TYPER_USE_RICH=0`), and `mcp` is absent from the tool list so the escaped copy reaches no description; the merge `030fd32b` diffed against each parent carries only the other side's rows in `test_layering.py` and `ROADMAP.md`, the `cli` row matches `cli.py` (`land_cover_fractions` at module level, `mcp.server` inside `mcp()`), `_core` rebuilt; `d8a1729d`: set equality implies the count, dropping `land_cover` from the derived `mesh` schema turns it False, 36 for `mesh` and 83 in all measured, the eight names pinned; section 12's paragraph and the five pinned "81" match the commits; check_citations exit 0; MCP suites plus layering 126 passed; ruff clean. Noted: the ctest 1031/1031 figure is in @developer's handback, not in `030fd32b`'s message; master's C++ is unchanged by this branch. Idea applied: ROADMAP row 37 says 328 lines.
