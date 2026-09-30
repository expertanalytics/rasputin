# Harness U1: queue, don't wait (with the U0 probes)

Status: design (@architect, 2026-09-30); U0 run, results in §2.5. No code. Source:
`docs/research/unattended-mode.md` (merged in #122), cited below as *the
design*, including *Ola's rulings (2026-09-30)* at its end. This file does not
restate the design; it fixes what U1 builds, exactly enough for `@tester` to
write the red suite.

**Goal.** While Ola is away, no guard asks: every act a guard would ask about
is refused with a reason, recorded in a queue, and waits for him. Ola alone
enters the mode, from his own terminal, for a stated duration. Independent of
the mode, `guard_push` stops missing the ref, remote, config and API writes
the design's §6 row *Publishing through the API* found unguarded.

**Scope ruling.** Ola, 2026-09-30: the three extra hooks of §3.7
(`AskUserQuestion`, `PermissionRequest`, `ConfigChange`) stay in U1; the scope
is as written here.

**Order.** Ola ruled on 2026-09-30 that U1 comes first ("U1 was the one I
wanted first"). The design's §8 placed U1 after generic-harness step 3
(history out of hook docstrings), which has not landed. The cost of the
inversion is one merge for step 3: U1 edits hook code, and the only docstring
lines it touches are the ones U1 makes false (§3.9).

---

## 1. Prior art: legacy and literature

*Literature.* Tooling; no novelty is claimed. Four established patterns are
combined, and each differs in one named way:

- **A state with an expiry, read on every use**: sudo's credential cache
  (`timestamp_timeout`, sudoers(5)). sudo's state *grants*; in U1 this flag only
  *restricts* (asks become refusals), so a forged or stale flag fails safe.
  U3 adds the design's one loosening, behind U2's push gate.
- **Refuse with a reason the caller sees**: Kubernetes validating admission
  webhooks (`allowed: false` plus a status message returned to the client).
  A PreToolUse `deny` has the same shape (design §1).
- **Approval held for a human**: GitHub Actions environments with required
  reviewers, where the job *waits*. Here the act is *refused and recorded*
  instead, because a waiting subagent stalls the whole run (design §1).
- **Confirmation read from `/dev/tty`, not stdin**: ssh and sudo read
  passwords from the controlling terminal so that a pipe cannot answer. Same
  trick, same limit: a pty wrapper defeats it, so it is a tripwire (design §2).

*Legacy.* Nothing to carry across. The legacy tree has no harness code:

```
$ git grep -l -iE 'unattended|caffeinate|isatty|/dev/tty' legacy-archive -- legacy
(no output, exit 1)
```

---

## 2. U0: the three probes

The design's §8 row U0 names three behaviours the documentation does not
settle: (a) a `deny` inside a background subagent, (b) a terminal behind `!`,
(c) `.claude/current-task/` writes in a worktree-isolated session. All three
run in **one sitting with Ola**, before the red suite, because (a)'s outcome
can stop U1 (§2.4).

### 2.1 Probe files (no production code)

Committed on the U1 branch by `@architect`, under `docs/increments/h3-probes/`
(precedent: `docs/increments/15-probes/`). They touch no `.claude/settings*.json`
and no live hook:

| File | Behaviour |
|---|---|
| `u0_probe_hook.py` | Reads the event from stdin. Appends one JSON line to the log (below; created if absent): `t` (UTC ISO), `hook_event_name`, `tool_name`, `agent_id`, `agent_type`, `permission_mode`, `cwd`, and the first 200 characters of `tool_input.command` or `tool_input.file_path`. Then, **only** for `PreToolUse` + `Bash` whose command contains `U0-PROBE-DENY`, prints a PreToolUse `deny` whose reason is `U0-REASON-7F3A: probe refusal. Quote this line in your handback, then run: echo U0-AFTER`. Otherwise prints nothing. Never exits non-zero |
| `u0-probe-settings.json` | Wires `u0_probe_hook.py` by absolute path on `PreToolUse` (matcher `*`), `PermissionRequest` (no matcher) and `SubagentStart` (no matcher) |
| `u0_tty_probe.py` | Prints and logs (same log, `hook_event_name: "tty-probe"`): `sys.stdin.isatty()`; whether `os.open("/dev/tty", os.O_RDWR)` succeeds (else the errno name); if it opened, writes `U0: type y and Enter within 20 s` to the tty and logs what `select` + `read` returned within 20 s. Exit 0: `y` read; 1: opened, no `y`; 3: no terminal |

**The log** is `/Users/skavhaug/projects/rasputin_scratch/u0-probe/log.jsonl`,
an absolute path in both scripts, so every session and subagent writes the
same file whatever its cwd. It does not exist before the sitting.

**Checked without Ola** (2026-09-30, from the Bash tool): the settings file
parses as JSON; the hook, fed a subagent `PreToolUse` Bash event containing
`U0-PROBE-DENY`, prints the `deny` with `U0-REASON-7F3A` and logs the line
with `agent_type`; an event without the marker, and non-JSON stdin, print
nothing, log, and exit 0; `u0_tty_probe.py` with stdin from `/dev/null`
prints `"stdin_isatty": false, "dev_tty": "ENXIO"` and exits 3. The log was
deleted afterwards.

The log is data, read by the session that ran the sitting; it is not a
channel (`.claude/REQUIRED-READING.md`, *Data, scratch and temp folders are not a channel*).

### 2.2 The sitting (Ola, once)

The temporary hook is **loaded by Ola's own launch flag**, so no
`.claude/settings*.json` is edited and there is nothing to revert. Ola, in his
own terminal, in the permission mode he uses when he leaves:

1. From the main checkout:
   `cd /Users/skavhaug/projects/rasputin && claude -w u0-probe --settings /Users/skavhaug/projects/rasputin/.claude/worktrees/unattended-u1/docs/increments/h3-probes/u0-probe-settings.json`
   (and `--permission-mode <mode>` if that is not the default).
2. At that session's prompt, probe (b):
   `! python3 /Users/skavhaug/projects/rasputin/.claude/worktrees/unattended-u1/docs/increments/h3-probes/u0_tty_probe.py`,
   and type `y` if asked.
3. At the same prompt, one line: *"Run probe U0 (a)+(c) as written in
   `/Users/skavhaug/projects/rasputin/.claude/worktrees/unattended-u1/docs/increments/h3-unattended-u1.md` §2.3, then print
   /Users/skavhaug/projects/rasputin_scratch/u0-probe/log.jsonl."* Ola may leave now. Any permission prompt that appears is itself a
   result: he answers it, and notes which it was.
4. When the session has printed the log, Ola closes it and tells the U1
   session "U0 done". The U1 session reads the log and `@architect` writes
   §2.5 (results) into this file.

**Fallback**, only if the log holds no `PreToolUse` line at all (the flag did
not load the hook): the same hook goes into the probe worktree's
`.claude/settings.local.json` (gitignored). That is an edit to the permission
system, so it needs Ola's fresh yes, once to add and once to remove.

### 2.3 What the probe session runs, and what decides each probe

The probe session's main agent, in order:

1. **(c), main session:** Write-tool write of `.claude/current-task/u0-probe.md`
   (one line), then a Bash write `echo u0 >> .claude/current-task/u0-probe.md`.
2. **(a), background subagent:** spawn `@tester` **in the background** with
   this brief: "U0 probe. (1) Write `.claude/current-task/tester-u0.md`
   (one line) with the Write tool. (2) Run `echo U0-PROBE-DENY`. (3) Do what
   the refusal says. (4) Hand back the refusal text verbatim and whether step 3
   ran." Wait at most 10 minutes for the handback.
3. Delete both `.claude/current-task/` files, then print the log.

| Probe | Confirms the design | Refutes it |
|---|---|---|
| (a) deny in a background subagent | The handback quotes `U0-REASON-7F3A`, the log has an `echo U0-AFTER` line with the same `agent_id` as the refused line, and no prompt appeared | *Stall*: no handback in 10 minutes, or a prompt appeared for the refused call. *Lost reason*: the subagent went on but did not see the text. *Ended*: the handback came, but no `U0-AFTER` line |
| (a′) `agent_type` (the R-B probe, generic-harness §7 step 8) | The refused line carries `agent_type: "tester"` | field absent or empty |
| (b) terminal behind `!` | `isatty` true, `/dev/tty` opens, and the `y` was read | either fails (`ENXIO` expected, as for the Bash tool, design §1) |
| (b′) do `!` commands pass through PreToolUse? | — (informational) | a `PreToolUse` line naming `u0_tty_probe.py` means they do |
| (c) current-task writes | no `PermissionRequest` line for either path, and no prompt | a `PermissionRequest` line naming `.claude/current-task/` |

The log's `permission_mode` field records which mode the sitting ran in, where
the event carries it.

### 2.4 What U1 does if a probe refutes

- **(a) Stall** stops U1 before the red suite: a deny inside a subagent is
  the mechanism, so the design's §1 inference is wrong and goes back to
  `@architect`. **Lost reason** or **Ended**: U1 ships as specified. The hook
  writes the queue itself, so the refusal still reaches Ola, via the recap and
  `away.py --back`; the subagent half of the queue reason is then advice that
  cannot arrive, and the result is recorded in §2.5.
- **(a′) absent**: the queue's `agent_type` is `null` for subagents too; the
  outcome is recorded for R-B.
- **(b)** changes no code. `away.py` refuses without a terminal in every
  outcome. It decides one line of rule text: Ola may run `away.py` with `!`
  only if (b) confirms **and** (b′) shows `!` bypasses PreToolUse, because
  U1's guard denies running `away.py` through any hooked Bash call (§3.6).
  Otherwise the rule says "in a separate terminal".
- **(c) refuted**: U1's `PermissionRequest` hook turns that prompt into a
  refusal at night, which would break every handoff file. U1 then adds allow
  rules `Write(.claude/current-task/**)` and `Edit(.claude/current-task/**)`
  to its settings change (Ola's yes, §3.10).

**Does the observation settle (c)?** Only for the configuration it was made
in. This session is worktree-isolated (its Bash tool refuses git commands it
cannot confine to `.claude/worktrees/unattended-u1`); its main session wrote
`.claude/current-task/session.md` and `architect-u1.md` there, and this
subagent rewrote `architect-u1.md` with the Write tool on 2026-09-30, and no
prompt reached it. That shows no stall for a worktree under
`.claude/worktrees/`, in whatever permission mode these sessions run. It does
not show *why* (the `.claude/worktrees` exemption, or the auto-mode classifier
approving), and it does not name the mode. So (c) no longer blocks U1 and
costs nothing extra: it rides along in the sitting, whose log records
`permission_mode`.

### 2.5 Results

Ola ran the sitting on 2026-09-30, 15:04-15:07 UTC, in the probe worktree
`.claude/worktrees/u0-probe`, launched with the `--settings` flag of §2.2. The
flag loaded the hook (the log has `PreToolUse` lines), so the fallback of §2.2
was not needed. The events record `permission_mode: "auto"`: the mode Ola
leaves in. The facts below come from the main session's reading of the log; the
log itself is data in `../rasputin_scratch` and is not cited here.

| Probe | Outcome | Evidence in the log |
|---|---|---|
| (a) deny in a background subagent | **Confirmed** | Background `@tester` (one `agent_id` throughout): `PreToolUse` Bash `echo U0-PROBE-DENY`, then 2 s later `PreToolUse` Bash `echo U0-AFTER` with the same `agent_id`, then the handback 7 s after that. `U0-AFTER` appears nowhere except the refusal reason, so the reason reached the subagent and it acted on it. Ola saw no prompt for the refused call |
| (a′) `agent_type` | **Confirmed** | Every subagent line, `SubagentStart` included, carries `agent_type: "tester"`. Main-session lines have `agent_id` and `agent_type` null |
| (b) terminal behind `!` | **Refuted**, as expected | The tty probe logged `stdin_isatty: false`, `dev_tty: "ENXIO"`, exit 3 |
| (b′) `!` through PreToolUse | **No**: `!` bypasses it | No `PreToolUse` line names `u0_tty_probe.py`. The one prompt Ola saw was a `PermissionRequest` for a later **Read** of the probe script by the probe session's agent (a path outside that worktree), not for the `!` command |
| (c) current-task writes | **Confirmed, in auto mode** | The main session's Write and Bash writes and `@tester`'s Write to `.claude/current-task/` produced no `PermissionRequest` line and no prompt |

**Consequences for U1** (per §2.4):

1. **(a) confirmed: U1 ships as specified.** A deny inside a background
   subagent neither stalls nor prompts, and its reason arrives, so the queue
   reason's instructions to a subagent (write the `ASK OLA:` line in the
   handback, do not retry) are read by the agent they address.
2. **(a′) confirmed:** the queue's `agent_type` and `agent_id` are populated
   for subagents and `null` for the main session, as §3.4 states. R-B
   (generic-harness step 8) can key on `agent_type`.
3. **(b) refuted: the rule text says "in a separate terminal"** (§3.12,
   unchanged). `away.py <duration>` run with `!` exits 3 without writing
   anything, which T12 already pins.
4. **(b′):** because `!` does not pass through PreToolUse, §3.6's guard never
   sees an `away.py` run with `!`; for entering, the terminal check of §3.3
   step 2 is the only barrier there, and (b) shows it holds. `--back` needs no
   terminal (§3.3), so Ola **can** end the mode with
   `! python3 tools/away.py --back` at a session prompt; §3.12 says so. No
   code changes.
5. **(c) confirmed in auto mode: no allow rules.** The allow rules §2.4
   names for a refuted (c) are not added to §3.11.

**The auto-mode caveat.** (c) is settled for auto mode only; other modes were
not run. U1 neither pins nor checks the mode:

- It cannot pin it: the permission mode is a property of each Claude Code
  session, set at launch, and `away.py` runs in Ola's terminal outside all of
  them.
- Checking buys nothing U1 needs. In a mode that prompts for writes, every
  file write prompts, not only those to `.claude/current-task/`, so
  unattended work in such a mode is not viable whatever U1 does; allow rules
  for one directory would not rescue it.
- It fails visibly, not silently: at night each such prompt becomes a
  `PermissionRequest` refusal (§3.7) with a queue line naming the tool and the
  path, which the recap and `--back` show.

If Ola starts leaving sessions in another mode, (c) is re-run in that mode
first (the §2.3 step 1 lines, with the probe settings), and the allow rules
of §2.4 are added to §3.11 if it refutes.

---

## 3. U1: the design

### 3.1 Components and data flow

```
Ola's terminal                        Claude Code (hooks, per tool call)
--------------                        -----------------------------------
tools/away.py <dur> --(/dev/tty y)--> <common>/harness/unattended.json  (flag)
tools/away.py --back ---------------> deletes flag; prints + archives queue
                                              |
                     tools/harness_mode.py  <-+-- read_mode()
                     (pure: decide(); IO: queue append)
                        ^            ^              ^
     guard_push.py -----+  guard_governance.py -+  guard_unattended.py (new)
     verdict: pass|ask|deny for its own patterns   AskUserQuestion, PermissionRequest,
                                                   ConfigChange: act only when mode != off
                        |
                        +--> <common>/harness/queue.jsonl  <-- tools/session_state.py (recap)
```

Boundaries. Each guard computes a **verdict** from the event alone and never
reads the mode. `harness_mode` maps (verdict, mode) to a decision and is
the only appender of `queue.jsonl`. `away.py` is the only writer of the flag
and of `windows.jsonl`, and `--back` archives the queue. `session_state.py`
only reads. No profile yet: the constants below live in
`harness_mode.py` until generic-harness step 4 moves them to
`.claude/profile.toml` (design §7).

**Locating the state dir.** Every caller passes its own checkout root, taken
from `Path(__file__)`: the hooks and `session_state.py` run from the live
checkout (`$CLAUDE_PROJECT_DIR`), and `away.py` from wherever Ola runs it;
all are in one repository. `state_dir(root)` runs
`git -C <root> rev-parse --path-format=absolute --git-common-dir` and appends
`harness`; it returns `None` if git fails. No environment variable redirects
it (tests copy the scripts into a temporary repository, §4).

Constants (rasputin values of the design's §7 profile keys): `OWNER = "Ola"`,
`MAX_HOURS = 72`, `BUFFER = 1.2`,
`KEEP_AWAKE = ("caffeinate", "-is", "-t", "{seconds}")`,
`FORGE_HOST = "github.com"`.

### 3.2 The flag

`<git-common-dir>/harness/unattended.json`. `away.py` creates the dir at mode
`0o700` and the file at `0o600`, written atomically (a temporary file in the
same dir, then `os.replace`). Content, per the design's §2:

```json
{"since": "2026-09-30T18:00:00+00:00", "until": "2026-10-01T03:36:00+00:00",
 "set_by": "away", "keep_awake_pid": 12345}
```

Times are `datetime.now(UTC).isoformat(timespec="seconds")`.
`read_mode(state: Path | None, now: datetime) -> Mode`, with
`Mode(state: Literal["off", "on", "broken"], until: datetime | None, detail: str)`:

| Condition, first match wins | Mode |
|---|---|
| `state` is `None` | `broken`, detail `git common dir not found` |
| flag absent | `off` |
| not a regular file; unreadable; not a JSON object; `since` or `until` missing, not a string, not `fromisoformat`-parseable, or without a UTC offset; `until <= since`; `until - since` > `MAX_HOURS × BUFFER` hours + 60 s | `broken`, detail naming the first failure |
| `now >= until` | `off` (expired; the file stays until `--back`) |
| otherwise | `on`, `until` set |

`keep_awake_pid` and `set_by` are informational and not validated.

### 3.3 `tools/away.py`

`python3 tools/away.py <duration>` and `python3 tools/away.py --back`. The
testable entry point is
`main(argv, *, root=None, now=None, open_tty=None, spawn=None, stop=None) -> int`;
the keyword seams default to the script's own checkout root,
`datetime.now(UTC)`, opening `/dev/tty`, `subprocess.Popen` and the
keep-awake stopper.

**Enter**, in this order; any failure writes nothing and starts nothing:

1. **Duration**: `^(?:(\d+)h)?(?:(\d+)m)?$`, at least one part, total
   minutes > 0 and ≤ `MAX_HOURS × 60` (the cap applies to the duration Ola
   types, before the buffer; ruling 3). Else exit 2 with the usage line.
2. **Terminal**: open `/dev/tty` read-write. On `OSError`, print to stderr
   `away.py needs your own terminal: /dev/tty is not available (<errno name>). Run it in a terminal window, not through the agent.`
   and exit 3. stdin is never read.
3. **Confirm** on the tty:
   `Unattended mode for <dur> x 1.2, until <local time> (<UTC ISO>). Guarded acts will be refused and queued.`
   plus, if a flag is on, `This replaces the window ending <until>.`, then
   `Type y to confirm: `. Only `y` or `yes` (case-insensitive, stripped)
   confirms. Else print `Not confirmed; nothing changed.` and exit 1.
4. **Battery**: if `pmset -g batt` succeeds and its output contains
   `Battery Power`, print
   `On battery: caffeinate -s holds only on AC power. Plug in, or the Mac may sleep.`
   A missing `pmset` is silent.
5. **Keep-awake**: stop the previous flag's `keep_awake_pid` (below), then
   start `KEEP_AWAKE` with `seconds = ceil(minutes × 60 × BUFFER)`,
   `start_new_session=True`, stdio to `DEVNULL`. If the executable is
   missing, `keep_awake_pid` is `null` and step 7 says
   `Keep-awake: not started (<why>)`.
6. **Write** the flag (§3.2) with `until = now + minutes × BUFFER`, then
   append `{"event": "enter", "since", "until"}` to `harness/windows.jsonl`
   (U2's SUSPECT mark reads it; design §4).
7. **Confirm on the terminal** (stdout), exit 0:
   `UNATTENDED until <local> (<UTC ISO>). Keep-awake: caffeinate pid <pid>. Guarded acts are refused and queued. End early: python3 tools/away.py --back`

**Stopping keep-awake**: only when the pid is an `int`, `ps -p <pid> -o comm=`
succeeds, and its basename is `caffeinate`; then `SIGTERM`. A reused pid is
never signalled.

**`--back`** needs no terminal: leaving loosens nothing (design §2), and the
guard still denies it to agents (§3.6). In order: read the flag in any
state; stop its keep-awake; delete it (absent is fine); append
`{"event": "back", "since", "until": <now>}` to `windows.jsonl` if a flag was
read; print the summary; archive the queue; exit 0. The summary:

```
Back. Unattended since <since> until <until as set> (ended <now>).
Queued while away (<n>):
  <branch>
    <at> <hook> <agent_type or main>: <act>
ASK OLA lines in .claude/current-task/:
  <file>: <line>
```

With no flag the first line is `No unattended flag was set.`; an empty or
absent queue prints `(none)`. Archive: `queue.jsonl` is renamed to
`queue-<UTC date>.jsonl`, or appended to that file and removed if it exists.

### 3.4 `tools/harness_mode.py`: `decide` and the queue

```python
Verdict = Literal["pass", "ask", "deny"]

@dataclass(frozen=True)
class Decision:
    kind: Literal["pass", "ask", "deny"]
    reason: str
    queue: bool          # the caller appends a queue entry

def decide(verdict: Verdict, reason: str, mode: Mode, act: str, why: str) -> Decision: ...  # pure
def queue_entry(event: dict, hook: str, act: str, why: str, now: datetime) -> dict: ...
def append_queue(state: Path | None, entry: dict) -> str | None: ...   # error text, or None
def guard(event: dict, hook: str, verdict: Verdict, reason: str, act: str, why: str) -> dict | None: ...
```

`guard` is the one call a PreToolUse guard makes: it resolves the state dir,
reads the mode, calls `decide`, appends to the queue when told to, and
returns the PreToolUse output (or `None` for pass).

| verdict | mode `off` | mode `on` | mode `broken` |
|---|---|---|---|
| `pass` | pass | pass | pass |
| `deny` | deny, guard's reason, no queue | same | same |
| `ask` | ask, guard's reason (today's text, unchanged) | **deny**, `QUEUE_REASON`, queue | **deny**, `BROKEN_REASON`, queue |

A subagent's `ask` stays `ask` while attended (R-B is not in U1).

- `QUEUE_REASON` is the design's §3 *queue reason*, verbatim, with `{until}`
  as `YYYY-MM-DD HH:MM UTC` and `{act}`, `{why}` from the guard.
- `BROKEN_REASON` replaces its first sentence with
  `Refused: the unattended flag cannot be read ({detail}), so every guarded act is refused until Ola runs python3 tools/away.py --back; {act} waits for Ola ({why}).`
  and keeps the rest.
- If `append_queue` returns an error, the sentence `The refusal is already
  recorded in the queue.` becomes `The refusal could not be recorded
  ({error}); your ASK OLA line is the only record.` The decision is still
  `deny`: a refusal never depends on the queue write.

**Queue entry**, one line in `harness/queue.jsonl`: the design's §4 fields
`at`, `branch`, `cwd`, `agent_type`, `hook`, `act`, plus `why` and
`agent_id` (the morning summary shows the reason, and tells two subagents of
one persona apart). `branch` is `git -C <event cwd> rev-parse --abbrev-ref HEAD`,
or `null`. `agent_type` and `agent_id` are `null` in the main session. `act`
is at most 1000 characters (the Bash command, or `<tool> <file_path>`), so
the line stays under `PIPE_BUF`. It is written with one `os.write` on an
`O_WRONLY | O_APPEND | O_CREAT` descriptor, mode `0o600`, creating the dir
if needed.

### 3.5 `guard_push`: the new patterns

Each new pattern yields verdict `ask`; `guard` then asks by day, and denies
and queues at night. The existing patterns are routed the same way, their
reasons unchanged. Matching is per **segment**: the command split on `&&`,
`||`, `;`, `|` and newlines, each segment tokenised with `shlex.split` (on
`ValueError`, `str.split`). A git pattern matches when a segment has a token
`git` and, after it, the subcommand token; options between them (`-C <path>`,
`-c <k=v>`) are skipped.

| Pattern | Asks for | Passes (read forms) | `why` |
|---|---|---|---|
| `git update-ref` | any | — | `update-ref moves a ref directly` |
| `git remote` | `add`, `set-url`, `rename`, `remove`, `rm`, `set-head`, `set-branches` | `-v`, `show`, `get-url`, no argument | `this changes where the remote points` |
| `git config` | any of `--add`, `--unset`, `--unset-all`, `--replace-all`, `--rename-section`, `--remove-section`, `--edit`, `-e`; or first positional `set`, `unset`, `rename-section`, `remove-section`, `edit`; or two or more positionals with none of `--get*`, `--list`, `-l`, `get`, `list` | one positional (`git config user.email`); `--get*`, `--list`, `-l`, `get`, `list` | `this writes git configuration (hooks path, remote URLs)` |
| `git symbolic-ref` | two or more positionals, or `-d` / `--delete` | one positional (`git symbolic-ref HEAD`, `--short HEAD`) | `symbolic-ref rewrites a symbolic ref` |
| `gh api` | `-X` / `--method` (as `-X POST`, `-XPOST`, `--method POST`, `--method=POST`) other than `GET`, case-insensitive; or, with no method given, any of `-f`, `-F`, `--field`, `--raw-field`, `--input` | no method and no field flag; explicit `GET` | `gh api with a writing method changes the forge` |
| `curl` to the forge | a token containing `FORGE_HOST`, and `-X` / `--request` other than `GET` or `HEAD`, or any of `-d`, `--data*`, `-F`, `--form`, `--json`, `-T`, `--upload-file` | a GET to the forge; any request to another host | `curl with a writing method to the forge` |

Positionals are tokens after the subcommand that do not start with `-`,
skipping the argument of `-f`/`--file`, `--blob`, `--type`, `--default`,
`--comment` (config) and `-m` (symbolic-ref). `gh api graphql` with `-f` is
asked: it is a POST, and a query cannot be told from a mutation by text.
The design's §6 row stands: this is a text scan, and a script file or
another HTTP client gets past it.

### 3.6 `guard_governance`

Three changes; the Bash arm's WRITES scan is otherwise unchanged.

1. **Routed through `guard`**: every hit is verdict `ask`, `why` =
   `a rule file changes`, `act` = `<tool> <path>` or the command. At night
   every rule-file edit is denied and queued, on every branch (the design's
   U1 row; U3 adds the allowed-and-logged arm).
2. **The self-protecting set** (design §3, ruling 4) joins the governed set,
   asked by day, denied at night: `GOVERNED` gains `tools/away.py`,
   `tools/harness_mode.py`, `tools/session_state.py`, `.claude/profile.toml`;
   `GOVERNED_PREFIXES` gains `.claude/skills/` and `.git/hooks/`; a new
   `GOVERNED_SUFFIXES = (".git/config",)` is matched by `endswith` only (the
   `tail ==` rule would govern every file named `config`). `.claude/hooks/`,
   the settings globs (which already match `~/.claude/settings.json`) and
   the rule files are governed today.
3. **Always denied, both modes, never queued** (design §2, *Agent writes the
   flag*), verdict `deny`:
   - Edit, Write or NotebookEdit whose path matches `(^|/)\.git/harness(/|$)`;
   - Bash naming `.git/harness` together with a WRITES match or any of
     `rm`, `touch`, `ln`, `mkdir`, `unlink`, `install` as a word;
   - Bash running `away.py`: a segment (as in §3.5) with a token whose
     basename is `away.py`, unless the segment's first token is one of
     `cat`, `less`, `head`, `tail`, `grep`, `rg`, `wc`, `diff`, `ls`, `git`,
     `ruff`, `mypy`, `pytest`, or `sed` without `-i`. A `script`, `expect`
     or `unbuffer` wrapper is thereby denied too.

   Reason: `Only Ola enters or leaves unattended mode, and only hooks and away.py write the harness state. Nothing is queued: this act is not an agent's to wait for.`

### 3.7 `guard_unattended.py` (new hook)

One script, dispatching on `hook_event_name`. When the mode is `off` it
prints nothing, for every event. Otherwise (`on` or `broken`) it appends a
queue entry and:

| Event | `act` | Output |
|---|---|---|
| `PreToolUse`, `tool_name == "AskUserQuestion"` | `AskUserQuestion: <first question, 200 chars>` | PreToolUse `deny` with the queue reason; `why` = `a question needs Ola; write it as the ASK OLA line` |
| `PermissionRequest` | `<tool_name>: <command or file_path, 200 chars>` | `{"hookSpecificOutput": {"hookEventName": "PermissionRequest", "decision": {"behavior": "deny", "message": <queue reason>}}}`; `why` = `a permission prompt needs Ola` |
| `ConfigChange` | `ConfigChange <source> <file_path>` | `{"decision": "block", "reason": "unattended: configuration changes wait for Ola"}` |

`@tester` checks the two output shapes against the `hooks` page (sections
*PermissionRequest* and *ConfigChange*) before pinning them; if the page
differs, the page wins, and this table is corrected in the same commit.

### 3.8 The recap (`tools/session_state.py`)

One line before `== recap ==` when the flag is not `off`:

- `on`: `UNATTENDED until <YYYY-MM-DD HH:MM UTC>. Guarded acts are refused and queued: record each as an ASK OLA line and continue.`
- `broken`: `UNATTENDED FLAG UNREADABLE (<detail>): every guarded act is refused. Ola: python3 tools/away.py --back.`
- expired flag still on disk: `Unattended mode ended at <until>. Ola: python3 tools/away.py --back prints and archives the queue.`

Inside the recap, after *Waiting on Ola*, whenever `queue.jsonl` holds a line
(in any mode: the queue outlives the window until `--back`):
`Queued while unattended (<n>):`, then the newest 10 as
`  <at> <branch> <hook> <agent_type or main>: <act, 100 chars>`, then
`  ... <n-10> older` if more, and `  (<k> unreadable lines)` if any.

Then, per the design's U1 row and §6 (*An uncommitted rule edit in the live
checkout*): `Uncommitted rule-file changes:` and `  <worktree> <path>` for
each path in `git status --porcelain` of each worktree in
`git worktree list --porcelain` for which `governed()` (imported by path from
`.claude/hooks/guard_governance.py`) is true; at most 10 lines, then
`  ... <n> more`. Omitted when there are none. Each added line is under 200
characters, so the recap stays far below the 10,000-character cap.

### 3.9 Failure direction

- A flag that cannot be read, or a state dir that cannot be found, is
  `broken`: every `ask` becomes a `deny` (design §3, *Failure direction*).
- A guard that raises after parsing the event emits `deny` with the
  exception's type and message, in both modes: a crash must not turn an ask
  into a pass. Unparseable stdin still prints nothing, as today.
- A hook that cannot import `harness_mode` denies whatever it would have
  asked, with the import error as the reason.
- Docstring lines U1 makes false are corrected in U1: `guard_push`'s
  "`ask`, not `deny`" and `guard_governance`'s "DECISION IS `ask`, NEVER
  `deny`" each gain the unattended exception. No other docstring text moves;
  generic-harness step 3 still owns the history.

### 3.10 What U1 does not include

- The acknowledgement gate: `harness_acks.py`, `harness_review.py`, the ack
  listing and SUSPECT mark, the prompt-ack arm, `REWORK` lines (**U2**).
- Night-time rule-file edits (allowed and logged on a branch), the pending
  journal, the "not in force" context, the `SubagentStart` hook (**U3**).
- R-B, a subagent's `ask` becoming `deny` by day (generic-harness step 8).
- The profile, `profile.toml` (generic-harness step 4).
- Signed acks (backlog, ruling 9).
- Ola's keep-awake memory note, which is his to edit (design §9, Q2).

### 3.11 Settings changes (Ola's fresh yes, once, for the U1 PR)

`guard_push.py` and `guard_governance.py` keep their wiring. Added to
`.claude/settings.json`, each running
`$CLAUDE_PROJECT_DIR/.claude/hooks/guard_unattended.py`:

1. `PreToolUse`, matcher `AskUserQuestion`.
2. `PermissionRequest`, no matcher.
3. `ConfigChange`, matcher `user_settings|project_settings|local_settings|skills`.

No allow rules for `.claude/current-task/`: probe (c) confirmed in auto mode
(§2.5).

Not a settings file: `pyproject.toml`'s `[tool.mypy] files` gains
`tools/harness_mode.py` and `tools/away.py`, so the strict gate covers them.

### 3.12 Rule text (rules only)

`.claude/REQUIRED-READING.md`, *The harness*: the `guard_push.py` list gains
the six patterns of §3.5, the hook list gains `guard_unattended.py`, and this
paragraph is added (the design's §7 text, cut to what U1 enforces):

> **Unattended mode.** Only Ola enters it, with `python3 tools/away.py
> <duration>` in a separate terminal, and it ends at the time the flag states
> or on `away.py --back` (which also works with `!` at a session prompt).
> No agent runs `away.py` or writes
> `<git-common-dir>/harness/`. While it is on, a guarded act is refused
> rather than asked, and the refusal is already queued: do not retry it or
> work around it; record it as an `ASK OLA:` line (main session: in
> `session.md`; subagent: in its handback) and continue with other work.
> When Ola says he is leaving, ask him how long, and ask him to run
> `away.py` with that duration.

"in a separate terminal" stands: probe (b) found no terminal behind `!`
(§2.5). The `--back` clause rests on probe (b′), `!` bypassing PreToolUse,
and on `--back` needing no terminal (§3.3). The file is governed, so the
edit asks when `@developer` makes it.

---

## 4. Tests `@tester` writes (red, before any of §3 exists)

**Fixture.** A temporary git repository (`git init`, one commit) into which
the test copies, at the same relative paths, `.claude/hooks/guard_push.py`,
`guard_governance.py`, `guard_unattended.py`, `tools/harness_mode.py`,
`tools/away.py` and `tools/session_state.py`. Hooks run as
`[sys.executable, <copy>]` with the event JSON on stdin and `cwd` = the
temporary repository; the state dir is then `<tmp>/.git/harness`, and the real
repository's is never touched. Flag helpers write it as `on`
(`until = now + 1 h`), `expired` (`since = now − 2 h`, `until = now − 1 s`) and
each `broken` row of §3.2 (not JSON, a directory, naive time,
`until <= since`, a window of `86.4 h + 61 s`). Pure functions are imported by
path, as `test_session_state.py` does, with `now` injected. Events carry
`hook_event_name`, `tool_name`, `tool_input`, `cwd`; subagent events add
`agent_id` and `agent_type`.

File names: `tests/python/test_harness_mode.py`, `test_guard_push.py`,
`test_guard_governance.py`, `test_guard_unattended.py`, `test_away.py`, and
additions to `test_session_state.py`. None is invariant-critical in the
README's sense: no mutation round.

| # | What | Pinned outcome |
|---|---|---|
| T1 | `read_mode`, each row of §3.2 | the mode, and a `detail` for each `broken` row; `now == until` is `off` |
| T2 | `decide`, the 3 × 3 table of §3.4 | kind, `queue` flag, and reason: day `ask` keeps the guard's reason; night reason contains `Refused: unattended mode is on until`, `ASK OLA:`, `handback`, `Do not retry`; broken reason contains `cannot be read` and `away.py --back` |
| T3 | queue write | one line per refusal with the fields of §3.4; two refusals give two lines; a state dir made read-only gives a `deny` whose reason contains `could not be recorded` |
| T4 | `guard_push`, each command below, **flag off, on, expired** | off and expired: `ask` rows ask, `pass` rows print nothing. On: `ask` rows deny with the queue reason and add one queue line (`hook` = `guard_push`, `act` = the command); `pass` rows print nothing and add none |
| T5 | `guard_push`, broken flag, `git push` | `deny`, broken reason, one queue line |
| T6 | `guard_governance`, Edit/Write of `CLAUDE.md` and of each new member (`tools/away.py`, `tools/harness_mode.py`, `tools/session_state.py`, `.claude/profile.toml`, `.claude/skills/x/SKILL.md`, `<tmp>/.git/hooks/pre-commit`, `<tmp>/.git/config`), flag off and on | off: `ask`; on: `deny` + queue line. `src/app/config` and `docs/config`: nothing, both modes |
| T7 | `guard_governance`, always-denied acts (§3.6 item 3), flag off **and** on | `deny`, no queue line: Write `<tmp>/.git/harness/unattended.json`; Bash `echo {} > .git/harness/unattended.json`, `rm .git/harness/unattended.json`, `python3 tools/away.py 8h`, `tools/away.py --back`, `script -q /dev/null python3 tools/away.py 8h`. Nothing: `cat .git/harness/queue.jsonl`, `cat tools/away.py`, `git diff tools/away.py`, `pytest tests/python/test_away.py` |
| T8 | subagent event (`agent_id`, `agent_type: "tester"`), `git push` | off: `ask`; on: `deny`, and the queue line carries both fields |
| T9 | crash and bad input | `tool_input` a JSON list: `deny`, both modes, both guards. Stdin not JSON: nothing |
| T10 | `tools/harness_mode.py` deleted from the copy, `git push` | `deny`, reason names the import error |
| T11 | `guard_unattended`, each event of §3.7, flag off, on, broken | off: nothing, no queue line. On and broken: the output shape of §3.7 and one queue line |
| T12 | `away.py` with **no terminal**: subprocess of the copy, `start_new_session=True`, stdin `DEVNULL`, argument `8h` | exit 3; stderr contains `needs your own terminal`; no flag, no `windows.jsonl`. In process, `open_tty` raising `OSError(ENXIO)`: `spawn` never called |
| T13 | `away.py` duration, in process | `8h`, `90m`, `1h30m`, `72h` accepted (480, 90, 90, 4320 min); `72h1m`, `0h`, `0m`, empty, `8`, `8x`, `-1h`, `1.5h`: exit 2, no flag |
| T14 | `away.py 8h`, fake tty answering `y` | exit 0; flag with `until − since` = 9 h 36 min, mode `0o600`; `spawn` called with `["caffeinate", "-is", "-t", "34560"]`; flag's `keep_awake_pid` = the fake's pid; one `enter` line in `windows.jsonl`; stdout contains `UNATTENDED until`. Answers `n` and empty: exit 1, no flag, no spawn |
| T15 | `away.py` edge cases | `spawn` raising `FileNotFoundError`: flag written, `keep_awake_pid` null, stdout contains `not started`. Re-entering while `on`: `stop` called with the old pid, prompt contains `replaces the window` |
| T16 | the stopper, real processes | a `sleep 30` child's pid is not signalled (alive afterwards); a non-`int` pid is ignored |
| T17 | `away.py --back` | with a flag, 3 queue lines on 2 branches, and an `ASK OLA:` line in `.claude/current-task/session.md`: stdout groups by branch and lists the line; flag gone; `queue.jsonl` gone and `queue-<date>.jsonl` holds 3 lines; one `back` line in `windows.jsonl`; `stop` called with the pid. A second `--back` appends to the same archive. No flag: `No unattended flag was set.`, exit 0 |
| T18 | recap | on, broken, expired: the first line of §3.8 in each case; off: no such line. 12 queue lines: 10 shown and `... 2 older`; one non-JSON line: `(1 unreadable lines)`; no queue: no section. `CLAUDE.md` modified in the repository and in a second worktree (`git worktree add`): both listed; a modified `notes.txt` is not |

**T4's commands.** Ask: `git update-ref refs/remotes/origin/master HEAD`;
`git -C /x update-ref -d refs/heads/y`; `git remote add up u`;
`git remote set-url origin u`; `git remote rename a b`; `git remote remove a`;
`git config core.hooksPath x`; `git config --unset remote.origin.url`;
`git config --global --add a.b c`; `git config set a.b c`;
`git symbolic-ref HEAD refs/heads/x`; `git symbolic-ref -d HEAD`;
`gh api -X PUT repos/o/r/pulls/1/merge`; `gh api --method=POST repos/o/r/pulls`;
`gh api -XDELETE repos/o/r/git/refs/heads/x`; `gh api repos/o/r/pulls -f title=t`;
`gh api graphql -f query=q`; `curl -X POST https://api.github.com/repos/o/r/pulls`;
`curl -d x https://api.github.com/x`; `echo ok && git update-ref a b`; and the
existing `git push`, `gh pr merge 1`, `git commit --amend`, `git rebase main`.
Pass: `git remote -v`; `git remote get-url origin`; `git config user.email`;
`git config --get remote.origin.url`; `git config --list`;
`git symbolic-ref --short HEAD`; `gh api repos/o/r/pulls`;
`gh api -X GET repos/o/r/pulls -f state=open`; `curl https://api.github.com/x`;
`curl -X POST https://example.org/x`; `git status`.

---

## 5. Estimate

Production lines as `CLAUDE.md` §2 counts them (tests excluded):

| File | Lines |
|---|---|
| `tools/harness_mode.py` (new) | ~110 |
| `tools/away.py` (new) | ~120 |
| `.claude/hooks/guard_push.py` (segmenter, six matchers, routing) | +~70 |
| `.claude/hooks/guard_governance.py` (members, always-deny, routing) | +~40 |
| `.claude/hooks/guard_unattended.py` (new) | ~45 |
| `tools/session_state.py` (header, queue, uncommitted scan) | +~40 |
| `.claude/settings.json`, `pyproject.toml` | +~16 |
| **Total** | **~440**, against the ceiling of 700 |

Above the design's ~260 because this file fixes what that figure left open:
the per-pattern parsing of §3.5, and the exact `--back` and recap formats.
If the green step runs over 600, the first cut is the uncommitted-rule-file
scan of §3.8 (~20 lines), which moves to U3.

---

## 6. What needs Ola

1. **The U0 sitting** (§2.2): done 2026-09-30; the flag loaded the hook, so
   no settings file was edited. Results in §2.5.
2. **The U1 settings change** (§3.11, items 1-3), one fresh yes when
   `@developer` reaches it.
3. **Leaving sessions in auto mode**, the only mode probe (c) covers (§2.5).
4. **Prompts during the green step**: `@developer` edits governed files
   (`.claude/hooks/*`, `tools/session_state.py`, the new `tools/away.py` and
   `tools/harness_mode.py` once governed, `.claude/REQUIRED-READING.md`), so
   the green step runs while Ola is at the keyboard.
5. The push and the PR, as always.
