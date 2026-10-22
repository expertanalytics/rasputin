# Increment 90: the MCP server picks up a changed sources file without a restart

**Status: built, approved** (red `7fb24be8`, green `38940a48`; 49 counted lines against about 40); designed by `@architect`, 2026-10-22, on master `94a157b2`,
branch `worktree-mcp-reload-sources`. Short path (`docs/increments/README.md`,
*Small changes*): about 40 counted lines, Python only, no refine or mesh code,
so no `@perf` run; `@reviewer` may send it to the full path. Number 90 because
89 is the highest number in `docs/increments/` over every local branch,
remote-tracking ref (both repositories) and worktree under
`.claude/worktrees/`, checked 2026-10-22 (`git ls-tree --name-only <ref>
docs/increments/`); 89 is the private package's.

## 1. What is asked, and why

Ola, 2026-10-22: "yes, make the server reload sources". That morning he
registered four new AR50 land-resource maps in `~/.config/rasputin/sources.toml`
and had to restart `rasputin mcp --http`, which voids every download link it
had handed out (the signing key lives only in memory, increment 69, 6.4).

## 2. What is stale today, measured

Each tool call runs a fresh child `python -m tin_engine <command>`
(`src_python/tin_engine/mcp/server.py@94a157b2:176`), and the child reads the
catalogue anew (`src_python/tin_engine/cli.py@94a157b2:1906`). So the
`sources` tool, and every source named by its id, already follow the file on
the next call, over both transports. Nothing else reads it.

One thing does not: over HTTP, the **read roots** (the folders a call may read
besides the caller's own) are computed once in `build_app`
(`src_python/tin_engine/mcp/http.py@94a157b2:225`), and they hold the `path`
of every local source. Probe (scratch, 2026-10-22, the test driver's
`rasputin_server` with the fake child): a sources file naming `a.tif`; after
start, `b.tif` added; a `mesh` call with `dem=[<path of b.tif>]` is refused
("outside your folder and the server's data sources"). Then `a.tif` removed
from the file; a call reading `a.tif` by path still runs. AR50 maps are
feature sources, which `sources` lists by path for `--features`, so this is
the case Ola hit.

Stdio has no read roots (`read_roots=None`), so stdio needs no change.

## 3. Prior art: legacy and literature

*Literature.* The pattern is nginx's configuration reload
(<https://nginx.org/en/docs/beginners_guide.html>, "Starting, Stopping, and
Reloading Configuration"): the new file is checked first; on failure "the
master process rolls back the changes and continues to work with the old
configuration"; old workers "continue to service current requests until all
such requests are serviced". This design keeps both rules: a refused file
keeps the last good roots, and a call already admitted keeps the roots it
started with. It differs in the trigger: nginx reloads on a signal; here the
file's content is compared at each call's start, because the person editing
the file is often not the one who started the server, and a check costs 17 µs
(section 6). Nothing new is claimed.

*Legacy.* Nothing. `git grep -n -i -E "reload|mtime|st_mtime|hot.?reload"
legacy-archive -- legacy` returned one line,
`legacy/rasputin/py2js.py:54` (a three.js `TextureLoader`), unrelated: the
legacy tree had no server.

## 4. Blueprint

```
build_app(app, settings, tokens, ...)                    [http.py]
  environ = dict(os.environ)                             copied once: the process's environment
  shared = SharedRoots(environ, out_root)                reads now; ValueError stops start, as today
  async scope_for(ctx) -> CallScope
      roots = await shared.current()
      CallScope(read_roots=(work, store, *roots), ...)   rest unchanged

SharedRoots.current() -> tuple[Path, ...]                [http.py]
  async with self._lock:                                 one reload at a time
      stamp = await to_thread(_stamp, candidates)        candidates = user_source_candidates(environ), fixed
      if stamp == self._stamp: return self.roots
      try:
          new = await to_thread(_load)                   _shared_roots(environ), then overlap(out_root, new)
      except ValueError as exc:
          LOG.warning("sources file refused; keeping the data sources read before: %s", exc)
          self._stamp = stamp                            judged once: a broken file logs once
          return self.roots
      self._stamp, self.roots = stamp, new               set after the load: a cancelled reload is retried
      LOG.info("read roots now %s", new)
      return new

call_tool / _download_link / _upload                     [server.py]
  scope = await scope_of(ctx)                            scope_for becomes async; stdio's is an async def
```

- **`SharedRoots`** (public in `tin_engine.mcp.http`, in `__all__`):
  `__init__(environ: Mapping[str, str], out_root: Path)`; attribute
  `roots: tuple[Path, ...]`, the last good set; `async current()`.
  `__init__` reads synchronously and raises what start raises today: a
  `ValueError` for a refused user file or for an overlap with `--out-root`.
- **`_stamp`**: for each candidate, its bytes; `None` when it is not a file;
  the `OSError`'s text when reading fails. The content, not the mtime, so a
  backup copied back with `cp -p` and an edit within one mtime tick are both
  seen, and a `touch` alone reloads nothing.
- **`_load`**: today's `_shared_roots(environ)` and `overlap(out_root, ...)`;
  `overlap` returns words, which `_load` raises as a `ValueError`.
- **`build_server`**: `scope_for: Callable[[ServerRequestContext[Any]],
  Awaitable[CallScope]] | None`; the default stdio scope is returned by an
  `async def`. The three `scope_of(ctx)` call sites gain `await`.
- **`HTTP_INSTRUCTIONS`**: "A change to the server's sources takes a restart."
  becomes "A change to the server's sources applies from your next call."

## 5. Rules

1. **A changed file applies from the next call.** The check runs when a call
   arrives (`scope_for`), for every tool, before the read check.
2. **A call already running keeps what it started with.** Its `CallScope` is
   a frozen value holding its own roots tuple; a reload replaces
   `SharedRoots.roots`, never a scope handed out. This includes a call
   waiting for a job slot: its second check (`src_python/tin_engine/mcp/server.py@94a157b2:541`) uses
   the same scope. Its child read the catalogue when it started.
3. **A refused file keeps the last good roots** and logs one warning line,
   naming the file and the reason (the catalogue's own words: file, line,
   column). The same refused content is not judged again, so a broken file
   logs once, not once per call. A `mesh`, `fetch` or catchment call is refused by its child, even one naming only paths, since each reads the catalogue (`src_python/tin_engine/cli.py@94a157b2:886`), in the CLI's words naming the file; that is the caller's report. The last good roots serve the server's own read check. Never all user sources dropped by
   a broken file.
4. **A new source that overlaps `--out-root` is refused the same way** (rule
   3): the overlap rule of increment 69, 5.2 is kept on every reload, since a
   root inside the output root would open other users' folders.
5. **A missing file is no user sources**, as the CLI reads it
   (`user_source_file` returns `None`): the roots become the cache and the
   plugins' local sources. Not a refusal; the server and its children agree.
6. **Removing a source** removes its path from the roots on the next call; a
   call reading it by path is then refused, by id refused by the child as
   unknown. Stored outputs and links are untouched: a link names a file in
   the user's `store` (a copy under its sha256), which is always a read root,
   and the signing key is never reloaded. Links survive a reload.
7. **Tokens stay read once**, at start, from `RASPUTIN_MCP_TOKENS`; the
   environment of a running process does not change, and nothing here reads
   it again. The environment copy in `build_app` fixes `RASPUTIN_DATA` and
   the XDG folders for the process's life, as today.
8. **Plugins' files are re-read with each reload but do not trigger one**:
   they change only with a package install, which needs a restart for its code
   anyway.
9. **Thread safety under `--jobs`.** `--jobs` caps the children; every check
   and reload runs on the server's one event loop. The lock makes concurrent
   arrivals after a change do one reload, not several, and log one line; the
   file read and the catalogue parse run in a worker thread, so the loop never
   blocks on them. The roots are an immutable tuple swapped in one assignment.
10. **Start is unchanged**: a refused file or an overlap at start still stops
    the server before it serves, since there is no last good set.

## 6. Cost

Measured 2026-10-22 on Ola's `~/.config/rasputin/sources.toml` (8,561 bytes,
14 sources), the worktree's venv: the stamp (read and compare one candidate)
17 µs per call, mean of 1,000; a reload (`_shared_roots`, the plugins'
entry points included) 3.1 ms, mean of 20. A tool call takes seconds, so no
timing gate. The constants assume a sources file of tens of kilobytes; the
largest checked is that 8.6 KB file.

## 7. Tests (`@tester`, red first)

New file `tests/python/test_mcp_reload_sources.py`. Fixture for U1 to U6: a
temporary folder with two 4-byte stub rasters `datasets/a.tif` and
`datasets/b.tif`, a sources file `data/sources.toml` written by
`dataset_fixtures.write_sources` with local entries (635 bytes for both
entries, measured), and an explicit `environ` of `RASPUTIN_DATA=<tmp>/data`,
`XDG_CONFIG_HOME=<tmp>/cfg`, `HOME=<tmp>`, so the user's real file is never
read; `out_root=<tmp>/out`. Reload count is counted from the log: one
`"read roots now"` INFO record of logger `tin_engine.mcp.http` per reload
(`caplog`). A reload costs about 3 ms (section 6), so the six unit tests
together stay well under a second.

- **U1 added.** File with `a`; build; rewrite with `a` and `b`;
  `await current()` holds `b.tif` resolved. Red: `SharedRoots` does not exist.
- **U2 removed.** Then rewrite with `b` alone; `a.tif` is gone from the roots.
- **U3 broken, kept, logged once, healed.** Append `[source.broken` (the
  catalogue refuses it: "Expected ']' at the end of a table declaration (at
  line 17, column 15)" on the two-entry file, measured); two `current()`
  calls return the previous roots, one WARNING record naming the file; fixed
  again, the next call reloads.
- **U4 overlap kept.** Add a source whose path is inside `out_root`; the roots
  are unchanged and the warning holds "overlap".
- **U5 unchanged content reloads nothing.** Two calls with the file untouched,
  then `os.utime` on it: zero "read roots now" records after build.
- **U6 concurrent arrivals reload once.** After a change, ten `current()` in
  one task group: all return the new roots, one "read roots now" record.
- **U7 missing file.** Delete the file; the roots are the cache alone.
- **S1 over HTTP, the wiring** (the probe of section 2, made a test):
  `rasputin_server(..., data=<tmp>/data)` with `a`; add `b`; a `mesh` call with
  `dem=[<b.tif>]` runs the fake child once (`set_spec(count_file=...)`, as in
  `test_mcp_http_read_limit.py`); remove `a`; a call reading `a.tif` is refused
  before the child. Red today on both halves (probe). The existing over-HTTP
  tests take 0.05 to 0.16 s a call (measured, `--durations`).
- **S2 a link survives a reload.** In S1's server, a `download_link` made
  before the file changes still downloads (200) after a reload.
  Green at red (nothing reloads the key); kept as the pin of rule 6.

Changed: `tests/python/mcp_http_driver.py`'s docstring of `rasputin_server`
says `RASPUTIN_DATA` holds "while the app is built"; with the environment
copied at build, it holds for the server's life, which S1 relies on. No
existing test changes its figures. Not invariant-critical: no mutation round.

## 8. Documentation (applied on this branch with the design)

- `README.md`, "Serving over HTTP": "Changing the sources takes a restart."
  becomes the paragraph in the README now (a change applies from the next
  call; a running call keeps its sources; a file that cannot be read keeps the
  last good sources, logged, and a call naming a source by id is refused until
  it is fixed; links survive).
- `docs/increments/69-mcp-http.md`, 5.2: the read roots "computed once at
  start" now say they are read again when the sources file changes, with a
  pointer here.

Checked and unchanged: README "Tokens are read at start, and links die with a
restart" (still true: tokens and the key are not reloaded); increment 37,
"the schema is read from the CLI at start" (the commands, not the sources);
increment 72 (nothing on sources); increment 78, a link from before a restart
works as an input (unchanged).

## 9. Left out

- Watching the file (kqueue, inotify): a dependency or a thread for a check
  that costs 17 µs per call.
- A `reload` tool or a signal: the content check needs no one to remember it.
- Handing the children the last good catalogue while the file is broken: the
  child reads the file itself (rule 3), and pointing it elsewhere would break
  the file's relative paths (`src_python/tin_engine/sources.py@94a157b2:386`).

## 10. Size

`http.py`: `SharedRoots` with `_stamp` and `_load`, about 30 lines;
`scope_for` and `build_app`, about 4. `server.py`: the type, the stdio
`async def` and three `await`s, about 5. Instructions: 0 net. About 40
counted lines (`python3 tools/count_loc.py <base> <head>`), worst case 55.

## Review

- **Design review, round 1, 2026-10-22.** On `721dfbb89be10b359425543064cf53de78514c8b` (94a157b2..721dfbb8, 0 counted lines), @reviewer: APPROVED WITH WORDINGS; two wordings (a broken file refuses every `mesh`, `fetch` and catchment call, not only calls by id, in §5 and the README; the stamp stored after the load so a cancelled reload is retried), applied in the commit that adds this line.
- **Code review, round 1, 2026-10-22.** On `38940a482bfd053a4f127babd7aedd51488a68f0` (94a157b2..38940a48, 49 counted lines; red 7fb24be8, green 38940a48; full suite rerun on a clean tree 9,573 passed), @reviewer: APPROVED WITH WORDINGS; two wordings (the status line, ROADMAP row 90), applied in the commit that adds this line.
