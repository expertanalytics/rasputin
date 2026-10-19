# Harness h26: no idle stop while away

Status: **designed**, 2026-10-19, on `worktree-h26-away-stop` off master
`38ff5e5f`. Short path (`docs/increments/README.md`): about 20 counted
lines, no C++ file, no refine or mesh code. Source: the retrospective
`docs/retrospectives/2026-10-19-d.md` (branch `worktree-retro-1019d`, PR
#302), section 2 (1 h 45 min idle in the afternoon away window) and its
proposal 1, which Ola said yes to on 2026-10-19 at about 20:1x UTC.

**The picture.** While Ola is away, nothing but a running agent, a
background command or a monitor wakes the main session once its turn ends.
So when the main session tries to end its turn with nothing running, a
`Stop` hook sends it back once, telling it to start queued work or write
down why there is none.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature* (read 2026-10-19, Claude Code 2.1.303 installed,
`claude --version`).

- Claude Code hooks reference, code.claude.com/docs/en/hooks, section
  *Stop*. The input carries, beside the common fields, `stop_hook_active`,
  `last_assistant_message`, `background_tasks` and `session_crons`.
  Verbatim: "The `stop_hook_active` field is `true` when Claude Code is
  already continuing as a result of a stop hook." "The `background_tasks`
  and `session_crons` arrays let hooks distinguish "session is done" from
  "session is paused waiting for background work to wake it back up". Both
  arrays are present when the task registry is reachable and are empty when
  nothing is in flight or scheduled." A `background_tasks` entry's `type`
  is a label "such as `shell`, `subagent`, `monitor`, `workflow`,
  `teammate`, `cloud session`, or `MCP task`". Decision: `decision:
  "block"` with a required `reason`. Loop guard: "after stop hooks have
  continued the turn eight times in a row, Claude Code overrides the next
  block and ends the turn. The count … resets each time Claude calls a
  tool." Stop "Runs when the main Claude Code agent has finished
  responding"; a subagent's end fires `SubagentStop` instead.
- Same page, failure: `onFailure` "has no effect on these hooks: Stop,
  SubagentStop …: exit code 2 on these events sends Claude back to keep
  working". Any other non-zero exit, a crash or invalid JSON is a
  "non-blocking error. The action goes ahead"; a hook that reaches its
  `timeout` "renders no decision". So the only way this hook traps the
  session is by exiting 2 or printing a block, and both are under its
  control.
- The `ralph-wiggum` plugin of the Claude Code repository blocks Stop to
  re-feed a prompt until a completion phrase or an iteration cap (secondary
  accounts, e.g. claudefa.st/blog/guide/mechanics/ralph-wiggum-technique;
  the plugin source was not read). Same mechanism; it loops on purpose,
  this hook refuses once.

The retrospective's open point, "`stop_hook_active` is from memory,
unchecked", is settled by the first quotation. Whether it stays `true`
after the continued turn has called a tool is not stated; section 4's live
check observes it, and section 3 shows either answer is safe.

*Legacy.* Nothing: the legacy tree is the old meshing product, with no
harness.

```
$ git grep -l -i -e 'stop_hook' -e 'unattended' -e 'away' legacy-archive -- legacy
(no output, exit 1)
```

*The ideas file* (`../rasputin_scratch/ideas/ideas.md`, read 2026-10-19):
no entry on idle stops; nothing taken up or overlapped.

## 2. The change

In `.claude/hooks/guard_unattended.py` (governed: its green step needs Ola
at the prompt), a fourth event, `Stop`, handled before `act_of` and never
queued (a stop is not an act waiting for Ola):

```python
def stop_output(event: dict) -> dict | None:
    """Refuse one idle stop of the main session while unattended mode is on."""
```

It returns `None` (the stop goes through) in each of these cases, checked in
this order, and the block below otherwise:

| Case | Why it passes |
|---|---|
| `event.get("stop_hook_active") is not False` | Already sent back once in this turn; a missing field also passes |
| `event.get("agent_id")` is set | A subagent; `Stop` should not fire there at all, this is a second fence |
| `"background_tasks" not in event` | The task registry was unreachable: running work cannot be ruled out |
| `background_tasks` or `session_crons` non-empty | Something will wake the session: any entry of any type or status counts |
| `harness_mode.read_mode(harness_mode.state_dir(harness_mode.ROOT), datetime.now(UTC))` is not `on` | Ola is present (`off`, expired), or the flag is `broken` |

This is the flag `harness_mode.settle` and `tools/away.py` read: the shared git common dir, so every worktree sees one flag. `stop_output` never calls `harness_mode.settle`, which would queue. The `harness_mode` import moves above the `act_of` check.

`broken` passes on purpose: this hook nudges rather than protects, a broken
flag already makes every other guard refuse, and Ola may well be at the
desk then; blocking would cost him one extra model turn per prompt.

The block, `{"decision": "block", "reason": STOP_REASON.format(until=...)}`,
with `until` as `harness_mode.decide` writes it (`%Y-%m-%d %H:%M UTC`):

> Unattended mode is on until {until} and nothing is running: no agent, no
> background command, no monitor, no scheduled wake-up. If this turn ends
> now, the session stays idle until Ola returns. Before ending it: write any
> open question as an `ASK OLA:` line in `.claude/current-task/session.md`
> and do not wait for the answer; start the next item that needs no ruling
> (the `QUEUE:` line), or watch work that will finish, for example
> `gh pr checks <pr> --watch` run in the background. If no such work exists,
> set the `NOW:` line to `NOW: idle: <why>` and end the turn; this hook does
> not refuse it a second time.

Failure: `main` already wraps `settle` in `except Exception` and prints
`refusal(name, …)`, which is `None` for any event but the three it knows,
so a crash on `Stop` prints nothing and exits 0: the stop goes through.
Unparseable stdin already exits 0 silently. The hook never exits 2. The
`git rev-parse` inside `state_dir` is the only external call; the settings
entry's timeout bounds it, and a timeout lets the stop through.

The module docstring's first line names four events. It says that Stop is refused once per turn only in mode `on` with nothing running. It also says a crash or a missing `harness_mode` lets a stop through, and a `broken` flag lets it through too. The test module's docstring says four hooks, not three (docstrings are not counted).

Estimate: about 20 counted lines (`stop_output` about 12, the reason
constant about 8 as a parenthesised string, one dispatch line in `settle`).

## 3. Why it cannot trap the session

Each refusal needs `stop_hook_active` to be `false`. If the field stays
`true` for the rest of the continued turn, there is one refusal per turn.
If it resets after a tool call, the session can be refused again only after
it has called a tool and still has nothing running, which is the case the
hook exists for; a stop with no tool call in between always passes, and
Claude Code's cap of eight continuations in a row overrides any block
after that. A crash, a timeout or a missing field lets the stop through.

## 4. Tests

In `tests/python/test_guard_unattended.py`, through `run_script` on a copied
repository as the existing tests do; the event is `{"hook_event_name":
"Stop", "cwd": …, "stop_hook_active": false, "last_assistant_message": "Done.",
"background_tasks": [], "session_crons": []}` with the field under test
changed. Fixture figures: the flag window is `harness_fixtures.flag_on`'s,
from now to now plus 1 h; `flag_on` returns `until`, which the test formats
rather than restating.

1. Mode `on`, idle event: stdout is `{"decision": "block", "reason": …}`,
   the reason contains the flag's `until` as `YYYY-MM-DD HH:MM UTC`, `ASK
   OLA:`, `QUEUE:` and `NOW: idle:`; exit 0; `queue_lines` is empty.
2. Mode `on`, each pass case of section 2's table (parametrised):
   `stop_hook_active: true`; `stop_hook_active` absent; `agent_id: "a1"`;
   `background_tasks` absent; one `background_tasks` entry for each of
   `shell`, `subagent` and `monitor` (`{"id": "t1", "type": …, "status":
   "running", "description": "x"}`); one `session_crons` entry. Each: exit
   0, empty stdout, empty queue.
3. Modes `off`, `expired` and `broken`, idle event: exit 0, empty stdout.
4. Mode `on`, `tools/harness_mode.py` deleted (as the existing
   missing-module test does): exit 0, empty stdout; the stop goes through.

The existing `test_the_hook_is_silent_while_attended` and the stdin tests
stay as they are.

*Live check* after Ola adds the settings entry, recorded in this file: in a
short window Ola opens (`away.py 15m`), the main session ends one turn with
nothing running, and the transcript shows one `Stop hook` refusal with the
text above; it then ends a second turn after one tool call, and the
transcript shows whether `stop_hook_active` stayed `true` (no second
refusal) or reset (a second one). A turn ended while a Monitor or a
background command runs shows no refusal. Without this check, a renamed or absent
`background_tasks` field would silently disable the hook (the third pass
case), so it is part of acceptance.

## 5. The settings entry (Ola's edit)

`.claude/settings.json`, a new top-level key beside `SubagentStop` (Stop
takes no matcher):

```json
"Stop": [
  {
    "hooks": [
      { "type": "command", "command": "$CLAUDE_PROJECT_DIR/.claude/hooks/guard_unattended.py", "timeout": 10 }
    ]
  }
]
```

No `onFailure`: it has no effect on `Stop`. The timeout of 10 s assumes one
`git rev-parse` on a local repository (well under a second on this
checkout); a slower run lets the stop through.

## 6. Left out

- A count of refusals on disk, for retrospectives to measure the hook's
  effect: the transcript shows each refusal already.
- Reading `session.md` to let a stop through once its `NOW:` line says
  `idle:`: `stop_hook_active` and the continuation cap make it unneeded.
- Refusing on a `broken` flag (section 2).

## Review
- **Design review, round 1, 2026-10-19.** On `f466d769`, 0 counted lines (estimate about 20), @reviewer: APPROVED WITH WORDINGS. Hook fields `stop_hook_active`, `agent_id`, `background_tasks`, `session_crons` checked against the Claude Code 2.1.303 hooks docs; no path traps the session; away flag read path and the docstring wordings applied by the main session.
