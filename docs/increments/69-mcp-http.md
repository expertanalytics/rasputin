# Increment 69: `rasputin mcp --http`, results returned as download links

**Status: designed** by `@architect`, 2026-10-18, on `cb713c87` (branch
`worktree-mcp-http`). Ola took the six defaults this design starts from
(section 1). Revised for `@reviewer`'s design review round 1 (changes
requested: outputs filtered, read limit independent of existence, re-check
after the job slot, `fetch --refresh` refused, quota figure). Next:
`@reviewer`'s second design review, then `@tester`'s red suite (section 9). Questions for Ola: section 12, each with a default the design
already follows.

## 1. What Ola asked

Ola, 2026-10-18: "I think it would be cool to get the MCP server setup dialed
in, so we communicate back using uris. It is some scaffolding, but will pay
off." The defaults he took for the first step:

- `rasputin mcp --http` serves MCP over Streamable HTTP, beside the stdio
  server of increment 37, in this repository, in `src_python/tin_engine/mcp/`.
- It runs on Ola's machine or one Expert Analytics server, behind TLS (a
  tunnel or a reverse proxy); no cloud account.
- Each user has one static bearer token, read from an environment variable;
  OAuth comes later. The first clients are Codex and Claude Code, which both
  send a configured bearer token.
- Results do not travel inside MCP messages. Each finished output is stored
  under its sha256 in a folder per user. The tool result gives a URI for it:
  a signed HTTP download link that the same process serves. The link is a
  bearer capability: it opens one file and it expires. Small text (the run
  record, `--stats`, legends) still comes back inline.
- Links are valid for 1 hour, stored outputs are kept 7 days, and each user
  has 20 GB.
- Resumable long runs (the MCP Tasks extension) are a later step.

Also required in this step: **a read limit**. Today the server limits where a
call writes (`_outside_root`, `src_python/tin_engine/mcp/server.py@cb713c87:263`)
but not what it reads, so any call can read any file on the machine. Over
HTTP, a call's inputs come only from the user's own folder or from the data
sources the server is configured with. Also required: path traversal and links
in user folders, token comparison in constant time, how the signing key is
handled, and what is logged.

## 2. Prior art: legacy and literature

*Literature.* The protocol is the Model Context Protocol, revision
**2026-07-28**, which `modelcontextprotocol.io/specification/versioning` names
as current (read 2026-10-18). The sections this design uses:

- *Transports, Streamable HTTP.* One MCP endpoint that accepts POST. Each
  request is answered with JSON or with an SSE stream scoped to that request,
  so progress notifications come before the result. Revision 2026-07-28
  removes protocol sessions and the GET stream. Earlier clients
  (2025-03-26 to 2025-11-25) still use `Mcp-Session-Id`. Security, verbatim:
  "Servers **MUST** validate the `Origin` header on all incoming
  connections"; "When running locally, servers **SHOULD** bind only to
  localhost (127.0.0.1)"; "Servers **SHOULD** implement proper authentication
  for all connections". Cancellation: "Closing the SSE response stream
  **MUST** be treated by the server as cancellation of that request." The
  spec also says servers "**SHOULD** include the `X-Accel-Buffering: no`
  header" on SSE.
- *Basic, Authorization.* "Authorization is **OPTIONAL**"; HTTP
  implementations "**SHOULD** conform", and conforming means OAuth 2.1 with
  Protected Resource Metadata (RFC 9728). Rules on token use: the token goes
  in the `Authorization: Bearer` header, "Access tokens **MUST NOT** be
  included in the URI query string", and "Invalid or expired tokens **MUST**
  receive a HTTP 401 response". **Departure:** a static token per user is
  not the OAuth flow, so the server publishes no resource metadata and no
  audience binding. That gives up the spec's guarantees on audience,
  rotation and revocation. Ola chose this for the first step; OAuth is the
  later step (section 8). The download link's signature does sit in a query
  string, but it is not an MCP access token: it opens one stored file until
  it expires, and it is never accepted at the MCP endpoint.
- *Server, Tools.* "A tool **MAY** return links to Resources", as
  `resource_link` content with `uri`, `name`, `mimeType` and `size`, and such
  links "are not guaranteed to appear in the results of a `resources/list`
  request". The download link is that kind of `uri`. `tools/list` "**MUST
  NOT** vary per-connection" but "**MAY** vary by the authorization"; here it
  is the same for every user. On handles, the non-normative *Stateful Tools*
  section says "For authenticated servers, a handle is a name, not a
  capability. The server should validate the caller's authorization against
  the handle on every call", and that retention "should be stated in the
  creation tool's description". Section 5.6 follows both. Its security list
  asks servers to "Rate limit tool invocations". **Departure:** this step
  has the existing jobs cap (calls queue) and the quota, but no rate limit,
  because the users are a handful of trusted colleagues (section 8).

*Capability URLs.* The download link follows W3C's "Good Practices for
Capability URLs" (First Public Working Draft, 18 February 2014): "Capability
URLs should be `https` URLs", "Capability URLs should expire", and an
expired one gets "either a `410 Gone` or a `404 Not Found` response". It
also warns that capability URLs "appear within application logs", which is
why section 6.5 turns the access log off. **Departure:** the draft asks that
a user can revoke a link. Here a link is revoked only by expiry or a server
restart, since the signing key lives only in memory (section 6.4). The link
is signed with HMAC-SHA256 (RFC 2104), the scheme of presigned object-store
URLs (for example S3's query-string authentication), applied here to files on
local disk. The store is content-addressed by sha256, as git stores objects.

*The SDK.* `mcp` **2.3.0**, the version `pyproject.toml` pins (`mcp>=2.3`)
and the newest on PyPI (`pip index versions mcp`, 2026-10-18). Facts this
design rests on, read from the installed package and checked by the probe
(section 3):

- `Server.streamable_http_app(*, streamable_http_path="/mcp", json_response,
  stateless_http, transport_security, host, auth, token_verifier,
  custom_starlette_routes, ...)` returns a Starlette app
  (`mcp/server/lowlevel/server.py`, line 732 of the installed package).
- `TokenVerifier` is a Protocol with `async verify_token(token) ->
  AccessToken | None` (`mcp/server/auth/provider.py`, line 124). The bearer
  middleware is installed only when `auth` (an `AuthSettings`) is given too.
  `AuthSettings` needs `issuer_url` and `resource_server_url`; the latter may
  be `None`, and then no metadata route is added.
- For an HTTP request, the handler's `ctx.request` is the Starlette
  `Request`, and `ctx.request.scope["user"]` is the `AuthenticatedUser`, whose
  `access_token.client_id` is what the verifier returned.
- The session manager records each session's owner and refuses another user's
  request on it (`streamable_http_manager.py`, lines 266-272).
- With `transport_security=None`, rebinding protection is switched on only
  for a loopback `host` (lines 753-758 of `server.py`); otherwise it is off.
- The request body limit is 4 MiB (`DEFAULT_MAX_REQUEST_BODY_SIZE`), the idle
  session timeout 30 minutes, at most 10,000 sessions.
- **The SDK client drops any SSE event over 1 MiB**
  (`DEFAULT_MAX_SSE_EVENT_SIZE = 1024 * 1024`,
  `mcp/client/streamable_http.py`, line 61). This is the message limit that
  rules out sending meshes inline (section 3).
- `X-Accel-Buffering: no` is set only on the 2026-07-28 code path
  (`_streamable_http_modern.py`). Responses to earlier clients do not carry it.

The SDK already depends on `starlette` (1.8.0 installed), `uvicorn` (0.54.0)
and `httpx2`. Starlette's `FileResponse` serves byte ranges
(`starlette/responses.py`, `accept-ranges: bytes`). The HTTP stack therefore
adds no new distribution, only two declarations (section 5.7). Signing and
token comparison use the standard library (`hmac`, `hashlib`, `secrets`);
`itsdangerous` is not installed and is not added. None of these is on the
prohibited list in CLAUDE.md §2.

*Clients.* Codex: in `config.toml`, a `[mcp_servers.<name>]` entry with `url`
and `bearer_token_env_var`, or `codex mcp add <name> --url <url>
--bearer-token-env-var <VAR>` (developers.openai.com/codex/mcp). Claude
Code: `claude mcp add --transport http <name> <url> --header
"Authorization: Bearer <token>"`.

*Novelty.* None is claimed. This is an adapter to a published protocol and
to well-known practice for capability URLs.

*Legacy.* Nothing carries across. The archived tree has no HTTP server, no
signing and no auth:

    $ git grep -l -i -E "uvicorn|starlette|hmac|bearer|signed.?url|streamable|http\.server|flask|aiohttp" legacy-archive -- legacy
    legacy-archive:legacy/rasputin/web_visualize.py

The one hit is a printed hint to run `python -m http.server 8080` for the
three.js viewer (lines 131 and 178), not code that serves anything.

*In this tree.* Increment 37's server (`src_python/tin_engine/mcp/server.py`,
`schema.py`): one child process per call, the write guard, progress
heartbeats, `resource_link` blocks with `file://` URIs. `DerivedCommand.writes`
(`src_python/tin_engine/mcp/schema.py@cb713c87:128`) already finds every writable path argument,
defaults included. This step adds its counterpart for reads. The catalogue
(`sources.catalogue`, `sources.user_source_file`, `sources.plugin_source_files`)
lists the local datasets whose `path` the read limit must allow.
`resources.parse_bytes` parses `20G`.

## 3. Measured (`docs/increments/69-probes/`)

**The SDK's HTTP auth, per user** (`http_auth_probe.py`, run with the main
checkout's venv, `mcp` 2.3.0, uvicorn on a free loopback port). Two clients,
one per token, call a tool at the same time. The handler sees its own
caller's user through both `ctx.request.scope["user"]` and the SDK's context
variable:

    [('ola', 'req=Request user=ola cv=ola'), ('kari', 'req=Request user=kari cv=kari')]
    no token: 401 Bearer error="invalid_token", error_description="Authentication required"
    bad token: 401

**Mesh sizes on real runs** (Romsdalen meshes from `../rasputin_scratch`,
binary `.vtk` and `.gpkg` of the same mesh; the triangle counts are the VTK
`POLYGONS` line and `select count(*) from triangles`, which agree):

| mesh | triangles | `.vtk` bytes | bytes/triangle | `.gpkg` bytes | bytes/triangle |
| - | - | - | - | - | - |
| romsdalen_uniform_10m | 615,459 | 22,195,158 | 36.1 | 153,145,344 | 248.8 |
| romsdalen_slope_2m_25_35 | 3,096,202 | 111,542,331 | 36.0 | 779,407,360 | 251.7 |
| velhas_t5_lcon (land-cover shares on) | 1,571,955 | 82,901,806 | 52.7 | - | - |

So 1 MiB holds about 29,000 triangles of binary `.vtk` and about 4,200 of
GeoPackage, before base64 adds a third. Real meshes run from tens to
hundreds of MB, and links are the only option.

**Storing an output** (the 779,407,360-byte GeoPackage above, on battery,
APFS): sha256 by `hashlib.file_digest` 0.38 s, `shutil.copyfile` 0.18 s. A
mesh that size takes minutes to build, so storing it adds well under 1 %.

## 4. Blueprint: components and data flow

```
client (Codex, Claude Code)
  │ POST /mcp, Authorization: Bearer <token>
  ▼
Starlette app from Server.streamable_http_app            [SDK]
  ├─ TransportSecurity: Host and Origin allow-list from --public-url
  ├─ BearerAuthBackend → StaticTokens.verify_token        [http.py]
  │     → AccessToken(client_id=<user>, token=<user>) or None → 401
  └─ call_tool(ctx, params)                               [server.py, as today]
        scope = scope_for(ctx)        CallScope: cwd, write root, read roots, user
        validate → refused flags (fetch --refresh)?       [workspace.py]
                 → writes inside scope.write_root?        [workspace.py]
                 → reads inside scope.read_roots?         [workspace.py]
                 → quota (HTTP only)                      [store.py]
        take the job slot, then run the same checks again
        child: python -m tin_engine <argv>, cwd = scope.cwd,
               env without RASPUTIN_MCP_TOKENS
        outputs: stdout lines resolved against scope.cwd;
                 HTTP keeps only regular files inside work/
              → delivery.deliver(scope, outputs)        [Delivery protocol]
            stdio: file:// URIs (today's behaviour)
            http:  store.put → sha256 → links.sign → https URL
  GET /files/<user>/<sha256><suffix>?fn=&exp=&sig=        [http.py route]
        links.verify → store path → FileResponse (no auth header needed)
```

What each module does and what it may not do:

- **`schema.py`** stays pure. It gains `DerivedCommand.reads(arguments) ->
  list[tuple[str, str]]`: (argument, value) for every string the call gives
  to a parameter that is not a writable path. That means `str` and `path`
  parameters, each item of a repeatable one, and each string element of a
  tuple (`tolerance_near` is `(path, float)`). The list includes defaults
  for the same reason `writes` does.
- **`workspace.py`** (new) decides which paths are allowed; it starts no
  child and holds no state. `CallScope` is a frozen dataclass: `cwd: Path |
  None` (None means the child inherits the server's cwd, as stdio does
  today), `write_root: Path`, `read_roots: tuple[Path, ...] | None` (None
  means no read limit, which is stdio), `user: str | None`, and `refused:
  frozenset[tuple[str, str]]`, the (tool, argument) pairs that may not be
  set on this transport (empty on stdio; `{("fetch", "refresh")}` over HTTP,
  section 7). One function, `async refusal(scope, tool_name, arguments,
  writes, reads) -> str | None`, gives the first refusal in plain words, where
  `arguments` is the validated model's `model_dump()`. It absorbs
  `_outside_root`.
- **`server.py`** keeps the child, the progress and the shaping of results.
  `build_server` gains two injected parameters with stdio defaults:
  `scope_for: Callable[[ServerRequestContext], CallScope]` and `delivery:
  Delivery`. `Delivery` is a Protocol with one method, `async deliver(scope,
  outputs: list[Path]) -> list[Delivered]`. Stdio passes a fixed scope and
  the file-URI delivery, so nothing outside the module changes, and
  increment 37's suite passes unchanged (T10). `call_tool` runs `refusal`
  (and, over HTTP, the quota check) before waiting for the job slot, so a
  refused call answers at once, and **again after taking the slot**, because
  the tree can change while a call waits. `_outputs` gains the scope:
  each stdout line is resolved as `(scope.cwd / line).resolve()`, not
  against the server's cwd, and when `scope.read_roots` is not None (HTTP)
  only a line that resolves to a regular file inside `scope.write_root` is
  an output (section 5.4). On stdio the rule stays today's (anything that
  exists).
- **`store.py`** (new): `Store(user_dir)` with `put(path) -> Stored(sha256,
  size, path)` (hash, then copy to `store/<sha256><suffix>` unless that file
  exists, whose mtime is then touched; run in a worker thread), `usage() ->
  int` (bytes under the user's folder), and `sweep(now, keep_seconds)`
  (removes regular files under `work/` and `store/` older than the keep
  time). Each is unit-testable on a temporary folder without a server.
- **`links.py`** (new, pure): `sign(key, user, name, filename, expires) ->
  str` and `valid(key, user, name, filename, expires, sig, now) -> bool`.
  Here `name` is `<sha256><suffix>` and the signature is
  HMAC-SHA256 over `user/name/filename/expires`, hex. Plus `url(public_url,
  ...)`. No I/O and no clock of its own: `now` is passed in.
- **`http.py`** (new, the only module that imports Starlette or uvicorn):
  `HttpSettings` (a frozen Pydantic model: host, port, public URL, out root,
  quota), `tokens_from_env(environ) -> dict[str, bytes]` (user to sha256 of
  the token), `StaticTokens` (the `TokenVerifier`), the HTTP `scope_for`
  and `Delivery`, the `/files` route, `build_app(...) -> Starlette`, and
  `serve_http(...)`, which runs uvicorn.
- **`cli.py`**: `mcp` gains `--http`, `--host`, `--port`, `--public-url` and
  `--quota`. Without `--http`, nothing about the command changes.

Data and execution stay apart: `HttpSettings` and `CallScope` are values,
and the modules that act on them take them as arguments. The one piece of
hidden state is the signing key, made once in `build_app` and held in a
closure. Everything is async; the hashing and copying go to
`anyio.to_thread`.

## 5. Decisions

### 5.1 The user folder

A user's folder is `<out-root>/<user>/` (`--out-root`, default
`~/rasputin-output`, as today). Inside it:

- `work/` is the child's cwd and the only place a call may write. A
  relative `out` such as `mesh.gpkg` lands here, and so does
  `sources`' `out_root`, since the child gets `RASPUTIN_MCP_OUT_ROOT=<work>`.
- `store/` holds `<sha256><suffix>` files and is written only by the server.
  Because no call can write there, a call cannot change a stored file under a
  link already handed out.

User names come from the token configuration only, never from a request.
They are checked when the server starts (`^[a-z][a-z0-9_-]{0,31}$`), so a
name can never be `..` or contain a `/`. Folders are made when a user's
first call arrives.

### 5.2 The read limit (HTTP only)

Every value from `reads()` is resolved as the child would read it:
`(scope.cwd / Path(value).expanduser()).resolve()`, following `..` and every
link. The call is refused, with the argument, the value as given and the
place it resolves to, when that place is outside every one of
`scope.read_roots` and either

- the value contains `/` or starts with `~`, **whether or not the file
  exists**; or
- the value has neither, and the resolved place exists. A bare name resolves
  inside `work/`, so this case is a link in `work/` pointing out, or `..`.

So a bare value that names no file is passed through: `glo30`, `corine`,
`EPSG:3006`, a layer name. The CLI then treats it as a name, or refuses it
as a missing file. Any relative value without `..` resolves inside `work/`
and is allowed whatever it holds (a URL such as `https://host/x` resolves to
`work/https:/host/x`). Because a value with `/` is judged without asking
whether it exists, `../kari/work/out.geojson` is refused before and after
Kari's run writes it, and the answer tells no caller which server files
exist. This rule needs no table of which parameters take paths. It holds for
plugin commands as well, and for the `str` options that take either a name
or a file (`--dem`, `--land-cover`, `--features-legend`, `palette`'s name).

The read roots, computed once at start (a change to the sources file takes a
restart, which the HTTP instructions say):

- the user's `work/`;
- the tile cache (`cli_kit`'s default: `$RASPUTIN_DATA/cache`, else the XDG
  cache);
- the `path` of every `kind = "local"` entry in the catalogue (the user's
  sources file and the plugins' files), each resolved.

The server refuses to start when a shared root (the cache or a catalogue
path) lies inside `--out-root` or contains it, since that root would open
other users' folders.

Stdio keeps no read limit (`read_roots=None`). Its user is on the same machine
and passes paths of their own today; limiting them is not part of this step.

### 5.3 The write guard over HTTP

The existing rule (`writes()` resolved inside the root, `..` and links
followed) is kept with `write_root = work/`, and relative values are now
resolved against `scope.cwd`, not the server's cwd. An explicit `cache`
argument to `fetch` is a write, so it must also lie inside `work/`. The
default cache is the server's and is shared (section 7).

### 5.4 What a result carries over HTTP

- `RunResult.outputs`: paths relative to `work/`, so a later call can take
  them as input (catchment's `.geojson` as mesh's `--domain`).
- A new optional `RunResult.downloads: list[Download] | None` (None on
  stdio), with one entry per output: `path` (as in `outputs`), `sha256`,
  `bytes`, `url`, `expires` (ISO 8601, UTC).
- One `resource_link` per output: `uri` = the signed URL, `name` = the file
  name, `mime_type` from `MIME_TYPES`, and `size`.
- `record`, `stats`, `stdout` and `stderr` as today, inline.

An output, over HTTP, is a stdout line that resolves against `scope.cwd` to a
**regular file inside the user's `work/`** (`lstat` after `resolve`, less the
sidecars). Nothing else is stored or linked: a line naming `/etc/hosts`, a
catalogue file, `../kari/work/x`, or a directory is left in `stdout` only.

### 5.5 The download route

`GET <public-url>/files/<user>/<name>?fn=<filename>&exp=<unix seconds>&sig=<hex>`.

- `user` must be a configured user; `name` must match
  `^[0-9a-f]{64}(\.[a-z0-9]{1,8})?$`; `fn` must match `^[A-Za-z0-9._-]{1,128}$`;
  `exp` must match `^[0-9]{1,12}$`. Otherwise **404**.
- The signature is checked with `hmac.compare_digest`. Bad: **403**.
- `exp` in the past: **410 Gone**.
- The path is built only from the checked parts, as `<out-root>/<user>/store/<name>`.
  It must be a regular file (`lstat`, so a link there is refused) or the
  answer is **404** (swept, or never stored).
- On success, `FileResponse` with `Content-Disposition: attachment;
  filename=<fn>`, the MIME type, byte ranges, `Cache-Control: private,
  no-store` and `Referrer-Policy: no-referrer`.

The route needs no bearer token. The link is the capability, so a client's
shell (`curl -o mesh.gpkg <url>`) can fetch it.

### 5.6 A fresh link: the `download_link` tool

Links last 1 hour but stored outputs last 7 days, and a restart voids every
link. So the HTTP server lists one tool beyond the commands, `download_link`,
with one argument `path`. The path must resolve inside the caller's `work/`
and be a regular file. It gets the same quota check as a command (refused at
or over the quota). The tool stores the file (so a re-made file gets a new
sha256) and returns one `Download` and one `resource_link`. The caller's user
is checked on each call, as the spec's *Stateful Tools* section asks. The
tool is not listed on stdio. If a command is named `download-link`, the
server refuses to start.

### 5.7 Dependencies

The `mcp` extra becomes `["mcp>=2.3", "anyio>=4.9", "mcp-types>=2.3",
"starlette>=0.49.1", "uvicorn>=0.31.1"]`. Both are already installed with the
SDK; they are declared because `http.py` imports them, and increment 48's
test fails on any undeclared import. `starlette>=0.49.1` is the release that
fixed a flaw in `FileResponse`'s Range parsing (CVE-2025-62727, per
starlette.dev/release-notes; range support itself came earlier);
`uvicorn>=0.31.1` is the SDK's own floor. Nothing goes into `dev` beyond what `mcp` brings.

### 5.8 Command line

```
rasputin mcp --http [--host 127.0.0.1] [--port 8765] [--public-url URL]
             [--quota 20G] [--out-root DIR] [--jobs N]
```

- `--host` defaults to `127.0.0.1`, following the spec's SHOULD.
- `--port 8765` is a chosen number.
- `--public-url` is the address clients use, through the tunnel or proxy;
  links are built on it, and its host and origin join the Host and Origin
  allow-lists along with loopback. It defaults to `http://<host>:<port>`.
  It is required when `--host` is not loopback, and it must be `https`
  unless its host is loopback.
- `--quota` takes the same forms as `--memory-budget`.

Links (3600 s) and retention (7 days) are constants in `http.py`
(section 10). `--http` without `RASPUTIN_MCP_TOKENS` is refused at start, so
the server is never served over HTTP without auth.

## 6. Security rules (each pinned by a test, section 9)

### 6.1 Tokens

- `RASPUTIN_MCP_TOKENS` holds `user=token` pairs separated by whitespace.
- At start, the server refuses the variable when it is missing or empty, or
  when it has a bad name, a repeated name, a repeated token, or a token
  shorter than 32 characters. 32 is a chosen floor:
  `secrets.token_urlsafe(32)` gives 43 characters, 256 bits. The refusal
  names the variable and the user, never a token.
- Only the sha256 of each token is kept. The `AccessToken` the verifier
  returns carries the user name in both `client_id` and `token`, so the token
  is not held in each request's scope.
- `verify_token` hashes the presented token and compares it, with
  `hmac.compare_digest`, against **every** stored digest, without stopping
  at a match. The time taken therefore does not depend on which user
  matched, and comparing fixed-length digests hides the token's length.
- The variable is removed from the child's environment.

### 6.2 Paths

- Sections 5.1 to 5.3, and the route's own checks (5.5).
- Every check resolves links, and `..` and `~` are expanded first.
- Sweep and usage walk with `os.walk(followlinks=False)` and remove only
  regular files.

The read rule does not depend on whether a file exists (5.2), so a file
written later, by any user, cannot turn a passed value into a read of
another user's folder: a value that resolves there is refused already. The
checks run again after the job slot is taken (section 4), so a link made
while a call waits is caught.

**Limit (assumption):** the second check and the child's open are still
separate moments. A file swapped for a link in that gap is not caught. A
user's calls write only in their own `work/` (5.3), so only that user's own
other call (with `--jobs` above 1) or `fetch` into the cache could do it, and
no command makes a link (`git grep -n symlink src_python` finds only
`cli_kit.py`'s refusal to follow one). Someone else with a shell on the
server is outside what this design defends against.

### 6.3 Host and Origin

The server always passes `TransportSecuritySettings(enable_dns_rebinding_protection=True,
allowed_hosts=[<public host[:port]>, "127.0.0.1:*", "localhost:*", "[::1]:*"],
allowed_origins=[<public origin>, "http://127.0.0.1:*", "http://localhost:*",
"http://[::1]:*"])`. It never relies on the SDK's loopback-only default.

### 6.4 The signing key

- 32 bytes from `secrets.token_bytes`, made once in `build_app`.
- Held in memory only: never written to disk, never put in the environment,
  never logged.
- A restart voids every link; `download_link` makes new ones.
- There is no key file to guard or rotate; the price is that links do not
  survive a restart.

### 6.5 Logging

- uvicorn runs with `access_log=False`, because its access line would hold
  the query string, and so the signature.
- The server logs to stderr through `logging`, as logger
  `tin_engine.mcp.http`. Per call, one line: UTC time, user, tool, exit
  code, seconds, output bytes. Per download, one line: UTC time, user, the
  first 12 hex of the sha256, status, bytes.
- Never logged: the `Authorization` header, a token, a signature, a full
  link, argument values.
- At start: the users' names, the public URL, and the read roots.

## 7. Degeneracy and edge policy

- **Over quota:** usage counts `work/` and `store/`, so an output is counted
  twice (once as written, once stored). When a user's folder holds at least
  `--quota` at a call's start, or at a `download_link` call, the call is
  refused, in words that give usage, quota and the 7-day
  sweep. A call that starts under the quota may finish over it. A plain
  refusal follows from the rule; no mid-run kill.
- **Sweep:** for the calling user, at server start and at each call's start.
  It removes regular files older than 7 days by mtime in `work/` and `store/`.
  Storing an output again touches its stored file. Links to a swept file
  answer 404.
- **Same output twice:** the second `put` finds `store/<sha>` and only
  touches it. Two runs that write the same bytes share one stored file.
- **Two calls by one user at once:** both pass the quota check and may both
  write. The jobs cap (default 1) serialises them anyway.
- **`fetch` over HTTP:** it writes the shared server cache. The cache is
  outside any user's quota, and a user can fill the disk with large fetches.
  Accepted for a handful of trusted users and named here; the cache is the
  operator's.
- **`fetch --refresh` over HTTP is refused** before the child starts:
  "fetch: --refresh discards the shared tile cache for every user, so it is
  refused over HTTP; give --cache inside your folder for a fresh copy".
  `--refresh` deletes the source's cache (`src_python/tin_engine/cli.py@cb713c87:1460`),
  which is every user's. `fetch-stations --refresh` is kept: its `--out-dir`
  is a write, so it lies inside the caller's `work/`.
- **Cancel:** the client closes the SSE stream (2026-07-28) or sends
  `notifications/cancelled` (earlier clients). The SDK cancels the handler,
  and `_run_child` kills the child, as on stdio. Nothing is stored.
- **A cancelled or failed call** stores nothing and returns no link.
  Files it left in `work/` are swept like any other.
- **An output the CLI names that is a directory** (a `station-catchments`
  `out_dir`): each regular file inside `work/` that the command printed is
  stored. A directory itself is never stored.
- **Server paths in `stdout`, `stderr` and `record`** are passed through as
  the CLI wrote them, so they name the server's folders. The users are
  colleagues, so this is accepted (section 8).

## 8. Left out, with reasons

- **OAuth** (resource metadata, audience binding): Ola's later step.
- **Tasks** (resumable long runs): Ola's later step. A long run is one SSE
  stream with heartbeats every 5 s, as on stdio.
- **Uploads:** no tool takes a file's bytes. Inputs are the server's data
  sources or earlier outputs (catchment, then mesh). The first gap a remote
  user will meet is "mesh my own polygon"; see question 1.
- **Rate limiting and per-user fairness:** beyond the jobs cap, none (trusted
  users).
- **Revoking one link:** only expiry or a restart.
- **Storing once:** each output is held twice, in `work/` and in `store/`.
  A hard link would hold it once, but a later run that rewrites
  `work/mesh.gpkg` in place would then change the bytes under a live link,
  and its sha256 would no longer match. A move would take the file away
  from later calls that read it (catchment's `.geojson` as mesh's
  `--domain`). Copy-on-write clones are not portable. The cost is half the
  quota (section 10).
- **Rewriting server paths** in stdout, stderr and record: not done.
- **`X-Accel-Buffering` for older clients:** the SDK sets it only for
  2026-07-28 clients. Behind nginx, turn buffering off in the proxy
  (`proxy_buffering off`). The README section says so.
- **Containers and multi-host deployment:** one process on one machine.
- **The read limit on stdio:** not changed.

## 9. Tests (`@tester`), red first

All of these go in new `tests/python/test_mcp_http_*.py` files, with
`pytest.importorskip("mcp")`. None needs the C++ kernel except T9. Calls go
through `mcp_fake_child.py`, increment 37's fake child, which writes the
files it is told to. Over HTTP, the server under test is `build_app(...)`
served by uvicorn on a free loopback port in a thread, as in the probe
(about 1.5 s to start; one server per module fixture). Where the transport is
not the point, the tests call the pure functions directly. **The
invariant-critical suite is T3 and T4** (the read limit and the link). Its
mutation targets are: `<` against `<=` on expiry, `compare_digest` replaced
by `==` in `links.valid` (killed by T4's wrapper test), the read check skipped
for tuple elements, `resolve()` dropped, and the `/`-or-`~` clause dropped
(so the check asks only whether the file exists; killed by T3's missing-file
case).

- **T1 tokens from the environment** (`tokens_from_env`, pure). It parses
  two users. It refuses a missing or empty variable, `Ola=…` (capital),
  `../x=…`, a repeated name, a repeated token, and a 31-character token; each
  refusal names the variable and the user, and no refusal contains the token
  text. A 32-character token is accepted.
- **T2 the verifier.** The right token gives `client_id` = the user. A wrong
  token, an empty one, and a correct token with one character added all give
  None. Over HTTP: no header gives 401, a wrong token 401, and a right token
  `tools/list` with the same tool names for both users, including
  `download_link`. With `hmac.compare_digest` patched to count its calls, one
  verify against two users calls it twice (no early exit).
- **T3 the read limit** (`workspace.refusal` on a temporary tree, then once
  over HTTP). Refused, each naming the argument: an absolute path outside the
  roots (`/etc/hosts`); `../other-user/work/x.geojson`; `~/x` where `HOME`
  is a temporary folder holding `x`; a link inside `work/` pointing outside;
  a directory outside given as `--dem`; `tolerance_near` = `["/etc/hosts",
  5]`; an item of repeatable `features` outside; **`../kari/work/out.geojson`
  when that file does not exist**, with the same refusal words (less the
  value) as when it does; `/nonexistent/x.tif`; `..` as a bare value. Allowed:
  a relative file in `work/`; a relative path naming no file
  (`sub/new.geojson`); a local catalogue source's path (a sources file written
  in the test); a bare value naming no file (`glo30`, `EPSG:3006`). Over HTTP,
  a refused call leaves the fake child uncalled. **After the slot:** with
  `--jobs 1`, a first call holds the slot (the fake child waits on a file the
  test creates); a second call with `--domain in.geojson`, a regular file in
  `work/`, passes the first check and waits; the test replaces `in.geojson`
  by a link to `/etc/hosts` and releases the first call; the second call is
  refused and the fake child is called once. `fetch` with `refresh=true` over
  HTTP is refused, naming `--refresh`, and the fake child is not called; on
  stdio it is not refused. On stdio (`read_roots=None`),
  `/etc/hosts` as `--domain` is not refused (today's behaviour kept).
- **T4 links** (`links`, pure, with `now` passed). A signed link checks out
  at `now = expires - 1`. At `now = expires` it is expired (the check is
  `now < expires`). It fails with any one of user, name, filename or expiry
  changed, and with another key. With `links.hmac.compare_digest` wrapped to
  record its calls (so `links` must call it through the `hmac` module), one
  `valid` call calls it exactly once, with the presented `sig` as one
  argument. Over HTTP: a stored file downloads with the
  stored bytes, `Content-Disposition` filename and MIME type. A `Range: bytes=0-9`
  request gets 206 and ten bytes. A changed `sig` gets 403, a past `exp`
  410, an `exp` of `1e9`, `-1` or 13 digits 404, a `name` of `../x` or
  `%2e%2e` 404, an unknown user 404. The `../x` request is sent raw with
  `http.client` (`/files/ola/../x?...`), because the SDK's `httpx2` client
  normalises it to `/files/x` before sending. A stored file
  removed by hand gets 404. A link replaced by a symlink in `store/` gets 404.
- **T5 a call over HTTP end to end** (fake child writing `a.gpkg` of 2 MiB,
  over the client's 1 MiB event limit). `outputs` is `["a.gpkg"]`.
  `downloads[0].sha256` equals the file's. The `resource_link` `uri` starts
  with the public URL and downloads the same bytes. The structured content is
  under 64 KiB. The file is at `<out-root>/<user>/work/a.gpkg` and
  `<out-root>/<user>/store/<sha>.gpkg`. A second identical call adds no
  file to `store/`. A fake child that also prints `/etc/hosts`,
  `../kari/work/x` (Kari's file made by the test) and a directory in `work/`
  gets exactly one download, `a.gpkg`; the three lines stay in `stdout`.
- **T6 two users apart.** Ola's and Kari's calls, made at the same time, write
  to their own `work/`. Kari's `download_link` on Ola's path is refused
  (outside her `work/`). Ola's session id sent with Kari's token is refused
  (the SDK's owner check, pinned so an SDK change is seen).
- **T7 quota and sweep** (`store.py` on a temporary folder, `now` passed).
  `usage` counts `work/` and `store/` and follows no link. `sweep` removes a
  file whose mtime is 7 days and 1 s old, keeps one 7 days less 1 s old, and
  removes no directory and no link target. Over HTTP with `--quota 1M`
  and 1 MiB already in `work/`, the next call is refused in words naming the
  quota, and the fake child is not called.
- **T8 start-up refusals and settings.** `--http` with no tokens exits 2. With
  `--host 0.0.0.0` and no `--public-url` it exits 2. With
  `--public-url http://example.org` it exits 2 (https needed off loopback).
  With `--public-url https://mcp.example.org` a request with
  `Host: mcp.example.org` passes, and `Origin: https://evil.example` gets
  403. The child's environment (the fake child dumps it) has no
  `RASPUTIN_MCP_TOKENS`. A plugin command named `download-link` stops start.
- **T9 a real `mesh` over HTTP** on the micro tile
  (`tests/python/geotiff_fixtures.py`, 3 × 4 cells, no codec). It returns
  exit 0, one `.vtk` in `downloads`, and the bytes from the link parse with
  `vtkread`.
- **T10 stdio unchanged.** Increment 37's MCP suites pass without edits
  (`pytest tests/python/test_mcp_*.py`).
- **T11 logging.** With `caplog` on `tin_engine.mcp.http` and uvicorn's
  loggers, one call and one download with a known token and signature: no
  record contains the token, the `sig` value or `?fn=`. One line names the
  user and the tool.
- **T12 progress over HTTP.** A fake child that prints three stderr lines
  0.1 s apart: the client's `progress_callback` sees three messages before
  the result. Cancelling the call (closing the client mid-run) kills the fake
  child within 2 s; the fake child writes its pid, and the test checks that
  the pid is gone.

## 10. Estimate, speed, acceptance

- **Counted lines**, by `python3 tools/count_loc.py`: `http.py` about 130,
  `store.py` about 45, `links.py` about 25,
  `server.py` net about 50 (the scope and delivery seam, the second check,
  the output filter, `download_link`, `downloads`), `workspace.py` about 60
  (refused flags, the `/`-or-`~` clause, the start-up root check),
  `schema.py` about 15, `cli.py` about 30, `pyproject.toml` 0 (not counted).
  **About 355; worst case 510.** That fits **one PR**, under
  the 700 of CLAUDE.md §2.
- **No `@perf` run:** no refine or mesh code is touched. Storing adds about
  0.6 s per 780 MB output (section 3).
- **Constants and their scale:**
  - Link life 3600 s, retention 7 days, quota 20 GiB (all Ola's). The quota
    assumes catchment-scale outputs. The largest checked is 779 MB, a
    3.1-million-triangle GeoPackage. Each output counts twice (`work/` and
    `store/`, section 8), so 20 GiB holds about 13 of them
    (21,474,836,480 / (2 × 779,407,360) = 13.8). A
    basin-scale GeoPackage can exceed the quota alone, hence `--quota`.
  - Token floor of 32 characters (chosen).
  - Port 8765 (chosen).
  - The 64 KiB ceiling on T5's structured content is a test bound against
    the client's 1 MiB event limit, checked at a 2 MiB output.
- **Docs in the same PR:** the README's MCP section gets a "Serving over
  HTTP" part: the two client commands from section 2, the token variable,
  the proxy note (section 8) and the 1 h / 7 d / 20 GB rules. The `mcp`
  command's help gains the HTTP options.

## 11. Ideas file

Read `../rasputin_scratch/ideas/ideas.md` on 2026-10-18. Entries this design
takes up or overlaps:

- "The MCP server as a service: results exchanged without a shared
  filesystem". This is its first step. It takes up Streamable HTTP, short-lived
  signed URLs, content addressing by sha256, and the write guard becoming a
  per-user workspace. It leaves out MCP resources read over the same
  connection (meshes are too big for the client's 1 MiB limit) and object
  storage.
- "Intent layer over rasputin". Overlaps only in serving MCP; not taken up.
- "Atomic mesh outputs". Overlaps: a `.vtk` killed mid-write is never stored,
  because only a call that exits 0 is stored. Left as its own entry.
- "Two small fixes from 37's code review" (the heartbeat task not awaited).
  A late progress tick over HTTP is a message after the result. Left out of
  this design and handed back.

## 12. Questions for Ola

1. **Uploads.** Should a remote agent be able to send its own small input
   files, such as a polygon to mesh? This step has no way to do that: inputs
   are the server's data and earlier results. *Default: not in this step.
   The next step adds a tool that takes a small text file (up to 1 MiB, the
   client's limit) into the user's folder.*
2. **Who may run `fetch` over HTTP.** It downloads public elevation data into
   the server's shared cache, outside any user's 20 GB. *Default: yes, every
   user may, but `fetch --refresh` is refused over HTTP, because it deletes
   that shared cache for every user.*

