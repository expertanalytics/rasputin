# Harness h28: no away window over an open prompt; the main session hears the window

Status: **designed**, 2026-10-20, on `worktree-h28-away-prompt` off master
`d3e55b0a`. Full path (`docs/increments/README.md`): about 75 counted lines,
over the short path's 50; no C++ file, no refine or mesh code. Every file it
changes is governed (`.claude/hooks/`, `tools/away.py`), so the red and green
steps need Ola at the prompt, by day.

Sources: Ola's yes of 2026-10-19 (`away.py` refuses to start while a
permission prompt is open, plus a Mac notification), and proposal 3 of the
retrospective `docs/retrospectives/2026-10-20.md` (branch
`worktree-retro-1020`): granted pushes waited 1 h 54 min to 3 h 07 min
because the main session did not notice the window had opened.

**The picture.** Three small parts. (A) The permission hook writes down
each prompt it lets open; `away.py` reads the session transcripts to see
whether each one was answered, and refuses to start while one is not.
(B) A Mac banner when a prompt has waited about six seconds, so Ola sees
it with the terminal hidden. (C) The first tool call of the main session
after a window opens gets a note: until when, and which branches are
granted.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature.* Claude Code hooks reference, code.claude.com/docs/en/hooks,
read 2026-10-20 (Claude Code 2.1.304 installed, `claude --version`).

- *PermissionRequest*: "Runs when Claude Code is about to ask you for
  permission to use a tool." Input: "`tool_name` and `tool_input` fields
  like PreToolUse hooks, but without `tool_use_id`." Common fields:
  `session_id`, `transcript_path`, `cwd`; `agent_id` and `agent_type`
  when the hook fires inside a subagent.
- *Notification*: the `permission_prompt` type fires when "Claude needs you
  to approve a tool use … and the prompt has waited about six seconds";
  "each keystroke defers it"; "If you or a PermissionRequest hook answer
  sooner, Claude Code doesn't run `permission_prompt`." Input: `message`,
  optional `title`, `notification_type`. "Exit code and stderr are
  ignored."
- *PermissionDenied* "only fires in auto mode: it doesn't run when you
  manually deny a permission dialog". **No event fires when a prompt is
  answered or dismissed**; the reference names none.
- *additionalContext*: returned in `hookSpecificOutput` of `PostToolUse`
  (the reference's own example); "Claude reads the reminder on the next
  model request, but it doesn't appear as a chat message".
- `transcript_path`: "written asynchronously and may lag the in-memory
  conversation". The reference does not say which transcript a hook inside
  a subagent is given.

**What the harness can see, then.** No hook, transcript field or state
file says "a prompt is open now". What can be seen: the moment a prompt
opens (PermissionRequest), and, in the transcript, whether that tool call
got a result. A rejected prompt writes one too: in
`~/.claude/projects/-Users-skavhaug-projects-rasputin/912df417-….jsonl`
(2026-10-05 08:06:02) a `user` entry holds a `tool_result` with
`is_error: true`, the `tool_use_id` of the call, and the text "The user
doesn't want to proceed with this tool use". So *open* is defined as: a
recorded PermissionRequest whose tool call has no `tool_result` in its
transcript. Checked on this session's live main transcript (508 MB): 2316
`tool_use` blocks, 2316 `tool_result`s, 0 unanswered; one scan took
0.85 s.

*Legacy.* Nothing: the legacy tree is the old meshing product, with no
harness.

```
$ git grep -l -i -e 'permission' -e 'osascript' -e 'unattended' legacy-archive -- legacy
(no output, exit 1)
```

*The ideas file* (`../rasputin_scratch/ideas/ideas.md`, read 2026-10-20):
no entry on prompts at the window's start, banners or a window note.
"Governed-file asks answered without a prompt in auto mode" touches prompts
but is a separate question; it is not taken up.

## 2. (A) A record of each prompt, and the check in `away.py`

**The record**, in `.claude/hooks/guard_unattended.py`, `settle`: when a
`PermissionRequest` gets `decision.kind == "ask"` (Ola present, so the
prompt opens), append one line to `<state>/prompts.jsonl` with
`harness_mode.append_queue(state, entry, "prompts.jsonl")`:

```python
{"at": <now, ISO, seconds>, "session_id": ..., "agent_id": ... or None,
 "agent_type": ... or None, "transcript_path": ..., "tool_name": ...,
 "input_sha": sha256(json.dumps(tool_input, sort_keys=True)).hexdigest(),
 "act": <act_of's act, first 200 characters>}
```

Nothing is recorded in mode `on` or `broken` (the prompt is refused and
never opens). A failed append is silent: the prompt still opens.

**The check**, in `tools/away.py`, a pure function over files:

```python
def open_prompts(state: Path, since: datetime | None) -> list[dict[str, Any]]:
    """Recorded prompts at or after `since` whose tool call has no result."""
```

- `since` is the `since` of the last `enter` line in `windows.jsonl`
  (`None` if there is none: every line counts). While a window is on, the
  hook records nothing, so lines before it were judged by that window's
  own start.
- The transcript of a line: if `agent_id` is set and
  `<transcript_path without .jsonl>/subagents/agent-<agent_id>.jsonl`
  exists, that file; else `transcript_path`. This works whichever path
  Claude Code hands a subagent's hook (section 5, live check).
- In the transcript, take the **last** `assistant` entry's `tool_use` block
  with the same `name` and the same `input_sha` (the same hash of its
  `input`). The prompt is open when that block's `id` has no `tool_result`
  in any `user` entry, and also when no matching block is found (the
  transcript may lag). A missing transcript file means the session is gone:
  not open. Read each transcript once, only the lines that contain
  `"tool_use"` or `"tool_result"`, and skip lines that fail to parse.

**The refusal.** In `enter`, before the terminal is opened, unless
`--despite-prompts` is given: if `open_prompts` returns any, print to
stderr and return 1, with nothing written and no keep-awake started:

> away.py: a permission prompt is open and would wait unanswered all
> night. {for each: session {session_id[:8]} ({agent_type or "main
> session"}), since {HH:MM UTC}: {act}}. Answer it, then run away.py
> again. If that session is closed, add --despite-prompts.

With `--despite-prompts` the same lines are printed as a warning, and
entering goes on. `--back` never checks. Only the prompt check reads
transcripts; `enter` has no other new step.

## 3. (B) The Mac banner

**What it is for:** a prompt that waits while Ola is present but looking
elsewhere. It fires when Claude Code sends the `permission_prompt`
notification, about six seconds after a prompt opened unanswered (the
3 h 13 min of 2026-10-18 and the 67 min of 2026-10-19, retrospectives
`2026-10-19-merges-290-292.md` §2.3 and `2026-10-19-merges-296-299.md`
§1.1). While away, prompts are refused before they open, so it never fires
then.

**How:** a new hook, `.claude/hooks/notify_prompt.py`, about 15 lines, no
new dependency:

```python
def notify(event: dict, run=subprocess.run, platform=sys.platform) -> list[str] | None:
    """The osascript argv for a permission_prompt notification, run once; else None."""
```

It acts only when `notification_type == "permission_prompt"`,
`platform == "darwin"` and `/usr/bin/osascript` exists. It runs

```
/usr/bin/osascript -e 'on run argv' -e 'display notification (item 1 of argv)
  with title "rasputin: a prompt is waiting" subtitle (item 2 of argv)'
  -e 'end run' <message> <worktree name>
```

with the message and the name of `cwd`'s last path part passed as arguments, so no
text is spliced into AppleScript; `timeout=5`, output discarded, any
`OSError` or `TimeoutExpired` swallowed, exit 0 always. The 5 s assumes one
local `osascript` call, measured at 0.12 s on this Mac (2026-10-20).
Exit 0 from `osascript` does not prove a banner showed: macOS lists such
banners under Script Editor, which may need notifications allowed in
System Settings once (Ola's check, section 5).

## 4. (C) The window note to the main session

In `guard_unattended.py`, a new event, `PostToolUse` (any tool), handled
before `act_of` like `Stop`:

```python
def window_note(event: dict, harness_mode) -> dict | None:
    """Tell the main session once per window that unattended mode is on."""
```

It returns `None` when `agent_id` is set, when the mode is not `on`, or
when `<state>/noticed.jsonl` already holds a line with this `session_id`
and the flag's `since`. Otherwise it appends that line
(`append_queue(state, {"session_id": ..., "since": ...}, "noticed.jsonl")`)
and returns

```python
{"hookSpecificOutput": {"hookEventName": "PostToolUse", "additionalContext": NOTE}}
```

> Unattended mode is on until {until}: Ola is away. Granted branches:
> {names, or "none"}. Each may be pushed, opened as a pull request and
> enqueued without asking, once `@reviewer`'s recorded verdict approves its
> head; do not wait for Ola's yes on these. Every other question is an
> `ASK OLA:` line.

`until` is formatted as `stop_output` does. `Mode` carries neither
`since` nor the grants, so the hook reads `<state>/unattended.json` once
more for both (`since`, and `push_grant.grants(flag)` for the names); a
flag that cannot be read then gives no note. A
crash prints nothing (`refusal` returns `None` for `PostToolUse`), so the
tool call is never affected. A new window (a new `since`) notes again.
Cost: one Python start and one `git rev-parse` per tool call, in every
session; mode `off` reads only the flag.

## 5. Settings entries (Ola's edit) and the live check

Two entries in `.claude/settings.json`:

```json
"Notification": [
  { "matcher": "permission_prompt",
    "hooks": [ { "type": "command", "command": "$CLAUDE_PROJECT_DIR/.claude/hooks/notify_prompt.py", "timeout": 10 } ] }
],
```

and under `PostToolUse`, beside the `Bash` entry, one with no matcher:

```json
{ "hooks": [ { "type": "command", "command": "$CLAUDE_PROJECT_DIR/.claude/hooks/guard_unattended.py", "timeout": 10 } ] }
```

The existing `PermissionRequest` entry already runs `guard_unattended.py`.

*Live check*, recorded in this file: (1) a prompt left open by a subagent
and one by the main session each make `away.py 15m` refuse, naming it;
answering it lets `away.py 15m` start (this also shows which transcript a
subagent's hook gets: `prompts.jsonl`'s `transcript_path`). (2) A banner
shows about six seconds after a prompt opens. (3) In that window, the main
session's next tool call shows the note once.

## 6. Tests

Fixtures are synthetic: transcripts of two to four JSON lines written in
`tmp_path`, `prompts.jsonl` of one to three lines, `windows.jsonl` of zero
or one `enter` line; the flag from `harness_fixtures.flag_on` (now to now
plus 1 h).

`tests/python/test_guard_unattended.py` (through `run_script`, as now):

1. Mode off, `PermissionRequest` for `Bash` `{"command": "git push"}`:
   stdout empty; `prompts.jsonl` has one line with the event's
   `session_id`, `transcript_path`, `agent_id`, `tool_name`, the
   `input_sha` of the input, and `act` `Bash: git push`.
2. Mode on, the same event: the refusal as today, no `prompts.jsonl`.
3. Mode on, `PostToolUse` of the main session: `additionalContext` contains
   `until` as `YYYY-MM-DD HH:MM UTC` and both granted names (flag with
   `grants: ["a", "b"]`); a second identical event: empty stdout; another
   `session_id`: the note; `agent_id` set: empty; mode off: empty; a flag
   with a new `since`: the note again.
4. Mode on, `tools/harness_mode.py` deleted, `PostToolUse`: exit 0, empty
   stdout.

`tests/python/test_away.py` (calling `main` with `root`, `now`, a fake
`open_tty` and `spawn`, as now), each with one recorded `Bash` prompt:

5. Transcript has the `tool_use`, no result: exit 1; stderr names the
   act and `--despite-prompts`; no flag, `open_tty` and `spawn` never
   called.
6. Result present (normal, and rejected with `is_error: true`): enters.
7. Transcript file missing: enters. No matching `tool_use`: refuses.
8. Two matching `tool_use`s, the first answered, the last not: refuses.
9. Prompt line older than the last `enter`'s `since`: enters.
10. `agent_id` set, the result only in `subagents/agent-<id>.jsonl`:
    enters; the result only in neither: refuses.
11. `--despite-prompts` with an open prompt: enters, and the warning names
    the act.

`tests/python/test_notify_prompt.py` (calling `notify` with a recording
`run`):

12. `permission_prompt`, message `say "hi"`, cwd `/x/h28-away-prompt`: one call; its
    argv ends with the message and `h28-away-prompt` unchanged, and no `-e`
    script contains the message.
13. `idle_prompt`, or `platform="linux"`: no call.
14. `run` raising `OSError` or `TimeoutExpired`: returns, no exception;
    the script on bad stdin exits 0.

Not the invariant-critical suite of any increment: no mutation round.

## 7. Size

About 75 counted lines (`python3 tools/count_loc.py <base> <head>`):
`open_prompts` and its transcript scan about 30, the refusal and the
`--despite-prompts` flag about 10 in `tools/away.py`; the prompt record
about 10 and `window_note` with its text about 18 in `guard_unattended.py`;
`notify_prompt.py` about 15 (counted: under `.claude/`, but Python). The
retrospective's new-module estimates ran 19 to 53 % low; the cap of 700 is
far off either way.

## 8. Left out

- A note when the window ends: Ola is then at the prompt.
- Clearing `prompts.jsonl` and `noticed.jsonl`: one line per prompt or per
  session and window, a few kilobytes a day; `open_prompts` reads only
  lines after the last window's start. Rotation can follow if they grow.
- Detecting a prompt by its six-second notification instead of
  `PermissionRequest`: it misses prompts answered within six seconds,
  which do no harm, but it also arrives late for the record; the
  transcript check is the same either way.
- A banner when `away.py` refuses: Ola is at that terminal.
