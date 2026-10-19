# Increment 72: MCP over HTTP, hardening after increment 69

**Status: built, approved** (red 7249a49f, green bd6b7054, test fix 93bf4cb8, 2 counted lines); designed by `@architect`, 2026-10-18, on `cef9f3f0` (master,
increment 69 merged as #290), branch `worktree-mcp-hardening`. **Short
path** (`docs/increments/README.md`, *Small changes*): about 12 counted
lines, no C++ file, no refine or mesh code.
Number 72 because 70 and 71 are held by designs on the worktrees
`lsm-ddd-routing` and `lsm-route-memory`.

## 1. What changes

Three follow-ups named by increment 69's code review
(`docs/increments/69-mcp-http.md`, *Review*, code review round 1), and one
documentation fix found while checking Ola's plan to run Codex against the
HTTP server (section 5).

### 1.1 Positional values go after `--`

`DerivedCommand.argv` (`src_python/tin_engine/mcp/schema.py@cef9f3f0:124`)
writes the command line in the click parameters' order, a positional value
as a bare word. A positional value that starts with `-` is then read by the
child as an option. Today the only positionals are `fetch`'s `source`,
`fetch-stations`' `source` and `palette`'s `name`, each a single required
`str`, so a value `--refresh` leaves the positional missing and the child
exits 2. But a plugin command with an optional positional, or two of them,
or one taking two values, would let a call slip a flag past the HTTP
refusals: `source="--refresh"` with `refresh=false` would run `fetch
--refresh` on the shared cache, and a two-value positional
`("--cache", "<tile cache>")` would turn a value the read limit allows (the
tile cache is a read root) into a write outside the caller's `work/`.

**The rule:** `argv` writes every option first, in the parameters' order,
then, when any positional is set, `--` and the positional values in the
parameters' order. Each item of a variadic positional (`nargs=-1`, a list)
is one word, as each element of a tuple already is. No positional set, no
`--`.

**Chosen over refusing read values that start with `-`**, because:

- `--` is the delimiter POSIX defines for exactly this (section 2), and
  click honours it (probe below). Every value keeps the meaning the caller
  gave it; nothing new is refused.
- A refusal would sit in `workspace.refusal`, which runs only where the
  read limit runs (HTTP). The argv is the same on stdio and HTTP, so the
  fix belongs where the argv is made, and it then holds on both.
- A refusal would also refuse a legitimate value, for instance a file in
  `work/` named `-old.geojson` given to a positional.

**Option values need nothing.** click takes the word after an option that
needs a value as that value even when it starts with `-`
(`click/parser.py`, `_get_value_from_state`, click 8.5.0), except for an
option whose value may be omitted (`_flag_needs_value`); no derived command
has one (probe below).

**The sidecars move into the model.** `_call` (`src_python/tin_engine/mcp/server.py@cef9f3f0:269`)
appends `--record <file>` after the argv and rewrites the word after the
first `--stats`. With positionals now last, an appended option would land
after `--` and be read as a positional, and `argv.index("--stats")` can find
an earlier option's value that happens to be `--stats`. Both are replaced by
building the argv from a copy of the arguments:
`tool.argv(arguments.model_copy(update={...}))` with `record_path` (when the
caller did not set it) and `stats` (when it is `-`). `model_copy(update=...)`
adds the updated keys to `model_fields_set` (checked, pydantic 2.13.5), so
`argv` writes them in their place. `mesh` has no positional today; this
keeps the rule true if it gains one.

**Probe** (`python -m tin_engine` from the main checkout's venv, click
8.5.0, cwd a scratch folder, 2026-10-18), what the CLI does with `--` for
each existing positional:

| command line | exit | what happened |
| - | - | - |
| `palette --out p5.toml -- corine` | 0 | `p5.toml` written |
| `palette --out p2.toml -- viridis` | 2 | "unknown legend viridis", as without `--` |
| `palette -- --help` | 2 | "unknown legend --help": a value, not the help |
| `palette --out p3.toml -- -x` | 2 | "unknown legend -x"; without `--`: "No such option: -x" |
| `fetch --bbox 0 0 1 1 --dry-run -- glo30` | 0 | the dry-run report |
| `fetch --dry-run -- corine` | 2 | "corine is not a catalogue source", as without `--` |
| `fetch -- --refresh` | 2 | "--refresh is not a catalogue source": a value, not the flag |
| `fetch-stations -- --refresh` | 2 | "Missing option '--out-dir'": `--refresh` taken as the source |
| `fetch-stations --out-dir o -- nosuch` | 2 | "nosuch is not a station list (nve-hrd)" |

So `--` changes no existing call's meaning, and after it a flag-like value
is a value. The parameter scan (`derive(app)` over the core app, no plugin
installed in that venv): the three positionals above, each `nargs=1`; no
parameter with `_flag_needs_value`.

### 1.2 A bad token pair is named by its position only

`tokens_from_env` (`src_python/tin_engine/mcp/http.py@cef9f3f0:109`) puts
the user part of a bad pair in its refusal. A pair written the wrong way
round, `token=user`, then prints the token at start-up. Checked, with a
32-character hex token whose first character is a letter (as `openssl rand
-hex 16` gives about 6 times in 16):

    RASPUTIN_MCP_TOKENS: the token of user <TOKEN> is under 32 characters

A `secrets.token_urlsafe` token, which the README recommends, mostly fails
the name check instead, and that refusal prints it too ("user <TOKEN> is
not a lowercase letter ...").

**The rule:** no refusal shows either side of a pair. Each names the
variable and the pair's position, counted from 1 in the order the pairs
are written, and the reason:

- `RASPUTIN_MCP_TOKENS: pair 2 is not user=token`
- `RASPUTIN_MCP_TOKENS: pair 2: the name is not a lowercase letter then up
  to 31 of a-z, 0-9, _ and - (a pair is user=token, the name first)`
- `RASPUTIN_MCP_TOKENS: pair 2 repeats the name or the token of an earlier
  pair`
- `RASPUTIN_MCP_TOKENS: pair 2: the token is under 32 characters`

The missing-or-empty refusal is unchanged. The docstring's "naming the
variable and the user" becomes "naming the variable and the pair's
position".

**Left as it is, named:** the start-up log line lists the user names
(`src_python/tin_engine/mcp/http.py@cef9f3f0:291`). A reversed pair whose "name" passes the name check
and whose "token" is 32 characters or more is accepted, and the log then
prints the token as a user name. That needs a real user name of at least 32
characters; names are capped at 32. Accepted.

### 1.3 T11 reads the SDK's server loggers too

T11 (`tests/python/test_mcp_http_calls.py`, from `is_the_servers` on)
reads records of `tin_engine` and `uvicorn` only. The `mcp` package logs
too. Measured (2026-10-18, a scratch copy of the test with `mcp` added to
both name lists, run against master): T11 passes, and the eight `mcp`
records of one call and one download all come from
`mcp.client.streamable_http`, the test's own client. The server's
`mcp.server.*` loggers wrote nothing on that path. They do write on
others: `transport_security.py` warns with the Host or Origin it refused,
and `streamable_http_manager.py` logs the first 64 characters of an unknown
session id. So adding the name alone would pass and measure nothing.

**The change (tests only, not counted):** `is_the_servers` takes
`mcp.server` and `mcp.shared` as well, never `mcp.client` (the test's
client, for the reason `httpx2` is left out). The test also sends one
request with a wrong token and one with `Origin: https://evil.example`, and
asserts that at least one record came from an `mcp.server` logger, so the
widened filter is shown to see something. The existing assertions (no
token, no `sig`, no key, no `?fn=`) then cover those records, and a new
one: the wrong token's text is in no record.

### 1.4 README: Codex's tool timeout

Codex stops waiting for a tool after `tool_timeout_sec`, default 60 seconds
(section 5). A large `mesh` takes minutes (increment 69, section 3). The
README's Codex snippet (*Serving over HTTP*) gains `tool_timeout_sec =
3600`, with one sentence, written in this design's commit: Codex's default of 60 seconds cuts off a long run. 3600 is a
chosen bound, one link life; the longest checked run is not measured here.

## 2. Prior art: legacy and literature

*Literature.* The delimiter is POSIX's: The Open Group Base Specifications
Issue 8, XBD section 12.2, Utility Syntax Guideline 10: "The first --
argument that is not an option-argument should be accepted as a delimiter
indicating the end of options", and the arguments after it are operands
even when they begin with `-`. The flaw it closes is CWE-88, "Improper
Neutralization of Argument Delimiters in a Command ('Argument Injection')"
(cwe.mitre.org/data/definitions/88.html, read 2026-10-18). The token rule
follows CWE-532, "Insertion of Sensitive Information into Log File": a
start-up refusal goes to stderr, which a service manager keeps as a log. No
departure from either. *Novelty:* none claimed.

*Legacy.* Nothing carries across; the archived tree builds no child command
lines and holds no tokens:

    $ git grep -l -E "subprocess|bearer|\"--\"|'--'" legacy-archive -- legacy
    (no files)

## 3. Tests (`@tester`), red first

All Python, none needs the C++ kernel. Figures are from the probes above.

- **R1 positionals after `--`** (`tests/python/test_mcp_schema.py`, the toy
  app). A toy command with a `--refresh` flag, one optional positional, one
  two-value positional and one variadic positional. For each case, the argv
  is parsed back by click itself (`command.make_context(name, argv[1:])`)
  and gives the values the call set: `{"first": "--refresh",
  "refresh": false}` parses to `first == "--refresh"` and `refresh` false;
  the two-value positional `["--cache", "x"]` parses to that tuple and no
  `cache` option; a variadic `["-a", "b"]` to those two items. The argv has
  options first, then `--`, then positionals; with no positional set it has
  no `--`. The eight rows of `TestArgvRoundTrip` whose toy command takes the
  positional `name` move it after `--` (for example
  `["make-thing", "--ascii", "--", "a"]`); the command with no parameters is
  still `["plain"]`.
- **R2 the real commands** (same file). For `fetch`, `fetch_stations` and
  `palette` from the core app, a call with the positional set to
  `--refresh` gives an argv that click parses back with the positional equal
  to `--refresh` and every flag at its default.
- **R3 sidecars before `--`** (`tests/python/test_mcp_server.py`, with the
  fake child, which already writes the argv it got). A toy `mesh` command
  given a positional: the argv the child gets has `--record` and `--stats`
  before `--`, each followed by its sidecar path; a real `mesh` call with
  `dem = ["--stats"]`, `out = "x.vtk"` and `stats = "-"` keeps that `dem`
  value and gets the stats sidecar after the real `--stats`. Today `argv` is
  `['mesh', '--out', 'x.vtk', '--dem', '--stats', '--stats', '-']` and
  `argv.index("--stats")` is 4, the `dem` value, which is rewritten (the red
  case; checked on master).
- **R4 token refusals by position** (`tests/python/test_mcp_http_tokens.py`,
  T1's parametrised refusal test, renamed to say "position"). Each refusal
  names the variable and `pair N` with the right N, and contains neither
  the name nor the token of the bad pair, checked as whole words (`ola` is
  a substring of ordinary words; pick names that are not, or match on word
  boundaries). New rows: a reversed pair `<32-hex token starting with a
  letter>=ola` (pair 1; neither side shown) and a reversed
  `token_urlsafe`-style pair `<token with a capital>=ola`; a bad second pair
  after a good first one gives `pair 2`.
- **R5 T11 widened** (section 1.3). Green on arrival: on master the
  bad-Origin request already gives `403` and one record,
  `mcp.server.transport_security` WARNING "Invalid Origin header:
  https://evil.example", and the wrong token `401` and none from the SDK
  (probe in a scratch copy of the test, 2026-10-18). It is a coverage
  change, not a red test, and the red commit's message says so.

Red on master: R1 (the flag-like and variadic cases, and the moved rows),
R2, R3 (the `--stats` value case), R4. **No mutation round**: nothing here
is the invariant-critical suite of an increment, and R1's round trip
through click is itself the check the review asked for.

## 4. Left out, with reasons

- **The heartbeat task not awaited** (`src_python/tin_engine/mcp/server.py@cef9f3f0:147`, `task.cancel()`
  with no await; ideas file, "Two small fixes from 37's code review"). Over
  HTTP a late tick would be a message after the result. Not one of 69's
  follow-ups; a deterministic red test for the race costs more than this
  increment. Handed back.
- **`--stats` swallows the next option** (ideas file). That is an option
  whose value may be omitted in the CLI; the MCP argv always gives `--stats`
  a value, so this increment does not touch it.
- **A download URL as an input.** Refused today (section 5); not changed.

## 5. Finding: Codex with no shell against `rasputin mcp --http`

Ola's plan: run Codex with no shell against the HTTP server, the exchange
limited to URIs. Nothing in increment 69 or this one stands in the way of
calls and results, with three things to know:

1. **Codex's tool timeout.** Codex's MCP settings
   (learn.chatgpt.com/docs/extend/mcp, read 2026-10-18): "`tool_timeout_sec`
   ... Default: `60`", "`startup_timeout_sec` ... Default: `10`". A mesh run
   longer than 60 s is cut off by Codex unless `tool_timeout_sec` is raised
   in its `config.toml`. Whether Codex resets that timer on progress
   notifications is not documented and was not checked. Section 1.4 puts the
   setting in the README; for tomorrow, set it by hand.
2. **Results are URIs; inputs are paths.** Each result carries its download
   links three ways: the `resource_link` blocks, `downloads[].url` in the
   structured content, and the same JSON as text
   (`src_python/tin_engine/mcp/server.py@cef9f3f0:313`), so a client that shows only text still
   sees them. A later call takes an earlier output by its `outputs` entry, a
   path relative to `work/` (`catchment.geojson` as mesh's `--domain`). A
   download URL given as an input is refused as "not a plain file path"
   (it holds `?`, `&` and `=`). With no shell, Codex cannot fetch a link
   itself; the links are for Ola's browser or `curl`.
3. **No uploads.** Codex cannot send its own polygon; inputs are the
   server's data sources and earlier outputs (increment 69, question 1).

Also checked: Codex sends no `Origin` header as a non-browser client, and
the SDK accepts a request without one (`mcp/server/transport_security.py`,
installed 2.3.0, line 75). Not checked: which MCP protocol revision Codex's
client speaks; the SDK serves earlier revisions with sessions (69, section
2), and T6 pins the session owner check.

## 6. Estimate

By `python3 tools/count_loc.py`: `schema.py` `argv` about +6 net (the
positional list, the `--`, the variadic items), `server.py` about 0 net (two
argv edits replaced by one `model_copy`), `http.py` about +2 net (the pair's
position; four messages reworded). **About 8; worst case 20.** One PR, short
path. No `@perf` run: no refine or mesh code. The README change is under
`README.md`, not counted.

**Constants:** `tool_timeout_sec = 3600` (chosen, a client setting, one
link life); no new constant in code.

## 7. Ideas file

Read `../rasputin_scratch/ideas/ideas.md` on 2026-10-18. Overlaps:

- "The MCP server as a service": this hardens its first step.
- "Two small fixes from 37's code review": the heartbeat await is left out
  (section 4) and handed back; the zero-bound fix is already in
  (`src_python/tin_engine/mcp/schema.py@cef9f3f0:66`, `is not None`).
- "`--stats` swallows the next option as its path": related (flag-like
  values), not taken up (section 4).

## Review

- **Combined round (short path), 2026-10-19.** Range `cef9f3f0..93bf4cb8`, 2 counted lines (estimate about 8), @reviewer: APPROVED WITH WORDINGS, applied by the main session. Full suite 8362 passed against this branch's code; 360 probed argv cases parse back exactly; user names stay in the log by design. Idea for later (on master too): two optional positionals where only the second is set are read as the first; refuse such a call.
